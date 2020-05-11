      subroutine w_mm_sort(ind_tau,w_rs,n1pbm,n2pbm,nd_nrr,n0_nrr,
     &                     n1pbmt,n2pbmt,nat1,nat2,isort,jsort,nrr_reds)
      use atom_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau,n1pbm,n2pbm,nd_nrr,n1pbmt,n2pbmt,
     &                       nat1,nat2,isort,jsort,n0_nrr,
     &                       nrr_reds(3,nqdiv*nat1*nat2)
      real*8, intent(out) :: w_rs(n1pbm,n2pbm,nd_nrr)
      integer :: k,k0,ind,iatom,ir0,ind_k,ir,jatom,ig,ind_ir0,j,i,
     &           jatom0,iatom0,l1,nl1,l,nl,je,ie,lm1,lm,km1,km,jat,jat0,
     &           iat,iat0,kmi,kma,kin
      real*8 :: pi2,v(3),phase
      complex*16 :: cc
      integer, allocatable :: km_row(:,:),km_col(:,:)
      real*8, allocatable :: ai(:,:,:),aj(:,:,:),a(:,:),b(:,:),c(:,:),
     &                       w0r(:,:,:,:,:)
      complex*16, allocatable :: w0(:,:,:,:,:),w1(:,:,:,:,:),
     &                           cphase(:,:,:,:),aa(:,:,:,:)
      pi2=pi+pi
      if(lsclp) then
        allocate(km_row(n_pbmt,ndim3_k(me_k+1)))
        allocate(km_col(n_pbmt,ndim3_k(me_k+1)))
        km_row=0
        km_col=0
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          do i=1,ncol_pw(k)
            j=refcol_pw(i,k)
            if(j<=n_pbmt) km_col(j,ind_k)=i
          enddo
          do i=1,nrow_pw(k)
            j=refrow_pw(i,k)
            if(j<=n_pbmt) km_row(j,ind_k)=i
          enddo
        enddo
      endif
      w_rs=0.d0
c ---- Store phase factors ---------------------------------------
      allocate(cphase(nqdiv,nat1,nat2,ngroup))
      do jatom0=1,nat2
        jatom=iat_sort(jatom0,jsort)
        do iatom0=1,nat1
          iatom=iat_sort(iatom0,isort)
          do ig=1,ngroup
            v=tshift(:,iatom,ig)-tshift(:,jatom,ig)
            do k=1,nqdiv
              phase=pi2*dot_product(pnt(:,k),v)
              cphase(k,iatom0,jatom0,ig)=dcmplx(cos(phase),sin(phase))
            enddo
          enddo
        enddo
      enddo
      if(ncmpl==1) then
        allocate(aa(2,0:maxpb,natom,ndim3_k(me_k+1)))
        call w_c_from_r_mm_prep(aa)
      endif
c -----------------------------------------------------------------
      do l1=0,lmpb(jsort)
        nl1=l1+l1+1
c ----- Store transformation matrix -------------------------------
        allocate(aj(nl1,nl1,ngroup))
        ind=l1*(nl1-2)*nl1/3
        do j=1,nl1
          do i=1,nl1
            ind=ind+1
            aj(i,j,1:ngroup)=u(ind,1:ngroup)
          enddo
        enddo
        do ig=1,ngroup
          do j=1,nl1
            do i=1,j-1
              phase=aj(i,j,ig)
              aj(i,j,ig)=aj(j,i,ig)
              aj(j,i,ig)=phase
            enddo
          enddo
        enddo
c ------------------------------------------------------------------
        do l=0,lmpb(isort)
          nl=l+l+1
c ----- Store transformation matrix -------------------------------
          allocate(ai(nl,nl,ngroup))
          ind=l*(nl-2)*nl/3
          do j=1,nl
            do i=1,nl
              ind=ind+1
              ai(i,j,1:ngroup)=u(ind,1:ngroup)
            enddo
          enddo
c ------------------------------------------------------------------
          if(ncmpl==1) allocate(w0r(nl,nl1,nat1,nat2,ndim3_k(me_k+1)))
          allocate(w0(nl,nl1,nat1,nat2,npnt))
          allocate(w1(nqdiv,nl,nl1,nat1,nat2))
          do je=1,ntle_pb(l1,jsort)
            do ie=1,ntle_pb(l,isort)
              if(ncmpl==1) w0r=0.d0
              if(ncmpl==2) w0=(0.d0,0.d0)
              do j=1,nl1
                lm1=l1*l1+j
                do i=1,nl
                  lm=l*l+i
                  if(lsclp) then
                    do ind_k=1,ndim3_k(me_k+1)
                      k=n3_mpi_k(me_k+1)+ind_k
                      do jatom0=1,nat2
                        jatom=iat_sort(jatom0,jsort)
                        km1=iopb(jatom)+ind_prod(je,lm1,jsort)-1
                        km1=km_col(km1,ind_k)
                        if(km1==0) cycle
                        do iatom0=1,nat1
                          iatom=iat_sort(iatom0,isort)
                          km=iopb(iatom)+ind_prod(ie,lm,isort)-1
                          km=km_row(km,ind_k)
                          if(km==0) cycle
                          if(ncmpl==1) then
                            w0r(i,j,iatom0,jatom0,ind_k)=
     &                        pw_sclp(km,km1,ind_tau,ind_k)
                          else
                            k0=2*km-1
                            call zcopy(1,pw_sclp(k0,km1,ind_tau,ind_k),
     &                                 1,cc,1)
                            w0(i,j,iatom0,jatom0,k)=cc
                          endif
                        enddo  !! over iatom0
                      enddo   !! over jatom0
                    enddo  !! over ind_k
                  else
                    do jatom0=1,nat2
                      jatom=iat_sort(jatom0,jsort)
                      km1=iopb(jatom)+ind_prod(je,lm1,jsort)-1
                      do iatom0=1,nat1
                        iatom=iat_sort(iatom0,isort)
                        km=iopb(iatom)+ind_prod(ie,lm,isort)-1
                        kmi=min(km,km1)
                        kma=max(km,km1)
                        kin=kma*(kma-1)/2+kmi
                        kin=kin-n_mpi_pbmt(me_b+1)
                        if(kin<=0) cycle
                        if(kin>ndim_pbmt(me_b+1)) cycle
                        if(ncmpl==1) then
                          do ind_k=1,ndim3_k(me_k+1)
                            w0r(i,j,iatom0,jatom0,ind_k)=
     &                        p_w_q_tau_mm(kin,ind_tau,ind_k)
                          enddo  !! over ind_k
                        else
                          k0=2*kin-1
                          do ind_k=1,ndim3_k(me_k+1)
                            k=n3_mpi_k(me_k+1)+ind_k
                            call zcopy(1,p_w_q_tau_mm(k0,ind_tau,ind_k),
     &                                 1,cc,1)
                            if(km<=km1) then
                              w0(i,j,iatom0,jatom0,k)=cc
                            else if(km>km1) then
                              w0(i,j,iatom0,jatom0,k)=conjg(cc)
                            endif
                          enddo  !! over ind_k
                        endif
                      enddo  !! over iatom0
                    enddo   !! over jatom0
                  endif
                enddo   !! over i
              enddo   !! over j
              if(ncmpl==1) then
                if(nproc_b/=1) call dgop(w0r,
     &                               ndim3_k(me_k+1)*nl*nl1*nat1*nat2,
     &                                   '  +',comm_b)
                w0r=w0r/2
                w0=(0.d0,0.d0)
                do ind_k=1,ndim3_k(me_k+1)
                  k=n3_mpi_k(me_k+1)+ind_k
                  call w_c_from_r_mm(w0(1,1,1,1,k),w0r(1,1,1,1,ind_k),
     &                               aa(1,0,1,ind_k),nl,nl1,nat1,nat2,
     &                               isort,jsort,l,l1)
                enddo
                if(nproc_k/=1) call dgop(w0,2*npnt*nl*nl1*nat1*nat2,
     &                                   '  +',comm_k)
              else
                if(nproc_k/=1) call dgop(w0,2*npnt*nl*nl1*nat1*nat2,
     &                                   '  +',comm_k)
                if(nproc_b/=1) call dgop(w0,2*npnt*nl*nl1*nat1*nat2,
     &                                   '  +',comm_b)
                w0=w0/2
              endif
c ------------------------------------------------------------------
              w1=(0.d0,0.d0)
              allocate(a(nl,nl1))
              allocate(b(nl,nl1))
              allocate(c(nl,nl1))
              do ind_k=1,ndim3_kk(me_k+1)
                k=n3_mpi_kk(me_k+1)+ind_k
                k0=i_kref(k)
                ig=k_group(k)
                call zone1_number(pnt(:,k),rb0,ndiv,ind)
                do jatom0=1,nat2
                  jatom=iat_sort(jatom0,jsort)
                  jat=ip(jatom,ig)
                  do j=1,nat2
                    if(jat==iat_sort(j,jsort)) jat0=j
                  enddo
                  do iatom0=1,nat1
                    iatom=iat_sort(iatom0,isort)
                    iat=ip(iatom,ig)
                    do i=1,nat1
                      if(iat==iat_sort(i,isort)) iat0=i
                    enddo
                    a=real(w0(:,:,iat0,jat0,k0))
                    b=matmul(ai(:,:,ig),a)
                    a=matmul(b,aj(:,:,ig))
                    c=imag(w0(:,:,iat0,jat0,k0))
                    b=matmul(ai(:,:,ig),c)
                    c=matmul(b,aj(:,:,ig))
                    w1(ind,:,:,iatom0,jatom0)=
     &                    cphase(k,iatom0,jatom0,ig)*dcmplx(a,c)
                  enddo  !! over iatom0
                enddo   !! over jatom0
              enddo   !! over ind_k
              deallocate(a,b,c)
c -----------------------------------------------------------------
              if(nproc_k/=1) call dgop(w1,2*nqdiv*nl*nl1*nat1*nat2,
     &                                 '  +',comm_k)
              do jatom0=1,nat2
                do iatom0=1,nat1
                  call fft3(ndiv(1),ndiv(2),ndiv(3),nl*nl1,
     &                      w1(1,1,1,iatom0,jatom0),1)
                enddo
              enddo
              w1=w1/dfloat(nqdiv)
              do ind_ir0=1,nd_nrr
                ir0=n0_nrr+ind_ir0
                ir=nrr_reds(1,ir0)
                iatom=nrr_reds(2,ir0)
                jatom=nrr_reds(3,ir0)
                do j=1,nl1
                  lm1=l1*l1+j
                  km1=ind_prod(je,lm1,jsort)
                  do i=1,nl
                    lm=l*l+i
                    km=ind_prod(ie,lm,isort)
                    w_rs(km,km1,ind_ir0)=w1(ir,i,j,iatom,jatom)
                  enddo
                enddo
              enddo  !! over ind_ir0
            enddo   !! over ie
          enddo   !! over je
          deallocate(ai,w0,w1)
          if(ncmpl==1) deallocate(w0r)
        enddo   !! over l
        deallocate(aj)
      enddo   !! over l1
      deallocate(cphase)
      if(lsclp) deallocate(km_row,km_col)
      if(ncmpl==1) deallocate(aa)
      end
