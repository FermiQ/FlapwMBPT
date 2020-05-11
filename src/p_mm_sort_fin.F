      subroutine p_mm_sort_fin(p_mm,nd_p2,n0_p2,n1pbm,n2pbm,nd_nrr,
     &                         n0_nrr,isort,jsort,nat1,nat2,nrr_reds,
     &                         r0_pairs,ig_pairs,ind_tau,pmm)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nd_p2,n0_p2,n1pbm,n2pbm,nd_nrr,n0_nrr,
     &                       isort,jsort,nat1,nat2,ind_tau,
     &                       nrr_reds(3,nqdiv*nat1*nat2),
     &                       r0_pairs(nqdiv,nat1,nat2),
     &                       ig_pairs(nqdiv,nat1,nat2)
      real*8, intent(in) :: p_mm(n1pbm,nd_p2,nd_nrr)
      real*8, intent(inout) :: pmm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                             ndim3_k(me_k+1))
      integer :: j,ig,iatom,jatom,ir,ir0,ij,l,lp1,ie,m,lm,iat,ind0,
     &           km,lm0,i,nn0,nn1,ir1,ind_r,ind_i,ind_j,ind_k,isize,jat,
     &           k,jnd0,km1,n0,n1,na
      real*8 :: y0(2*maxpb+1),t0(2*maxpb+1),phase,pi2
      integer, allocatable :: nd_p1(:),n0_p1(:)
      real*8, allocatable :: tmp(:,:),t(:,:),f0(:,:,:),f(:,:,:),
     &                       trf(:,:,:)
      complex*16, allocatable :: tc(:,:),p00(:,:,:)
      pi2=pi+pi
      allocate(nd_p1(nproc_b))
      allocate(n0_p1(nproc_b))
      call size_shift_par(n1pbm,nproc_b,nd_p1,n0_p1)
      nn0=n0_nrr
      nn1=nn0+nd_nrr
      allocate(trf(ndim3_kk(me_k+1),npnt,2))
      do ind_r=1,ndim3_kk(me_k+1)
        ir=n3_mpi_kk(me_k+1)+ind_r
        do k=1,npnt
          phase=pi2*dot_product(pnt(:,k),rr_3(:,ir))
          trf(ind_r,k,1)=cos(phase)
          trf(ind_r,k,2)=-sin(phase)
        enddo
      enddo
c ------------------------------------------------------------------
      if(ncmpl==1) then
        allocate(p00(n1pbm*nat1,n2pbm*nat2,ndim3_k(me_k+1)))
        p00=(0.d0,0.d0)
      endif
      n0=n_mpi_pbmt(me_b+1)
      n1=n0+ndim_pbmt(me_b+1)
      do jatom=1,nat2
        na=nat1
        if(isort==jsort) na=jatom
        do iatom=1,na
          allocate(f0(n1pbm,nd_p2,ndim3_kk(me_k+1)))
          allocate(t(n1pbm,nqdiv))
          do j=1,nd_p2
            t=0.d0
            do ir=1,nqdiv
              ir0=r0_pairs(ir,iatom,jatom)
              if(ir0>nn0.and.ir0<=nn1) then
                ir1=ir0-nn0
                t(:,ir)=p_mm(:,j,ir1)
              endif
            enddo
            if(nproc_k/=1) call dgop(t,n1pbm*nqdiv,'  +',comm_k)
            do ind_r=1,ndim3_kk(me_k+1)
              ir=n3_mpi_kk(me_k+1)+ind_r
              f0(:,j,ind_r)=t(:,ir)
            enddo
          enddo
          deallocate(t)
          allocate(tmp(n1pbm,n2pbm))
          allocate(f(nd_p1(me_b+1),n2pbm,ndim3_kk(me_k+1)))
          f=0.d0
          do ind_r=1,ndim3_kk(me_k+1)
            ir=n3_mpi_kk(me_k+1)+ind_r
            ig=ig_pairs(ir,iatom,jatom)
            tmp=0.d0
            do ind_j=1,nd_p2
              j=n0_p2+ind_j
              do l=0,lmpb(isort)
                lp1=l+1
                do ie=1,ntle_pb(l,isort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,isort)
                    lm0=l+m+1
                    t0(lm0)=f0(km,ind_j,ind_r)
                  enddo
                  call rotate1(t0,y0,l,u(1,ig),1)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,isort)
                    lm0=l+m+1
                    tmp(km,j)=y0(lm0)
                  enddo
                enddo
              enddo   !! over l
            enddo   !! over ind_j
            if(nproc_b/=1) call dgop(tmp,n1pbm*n2pbm,'  +',comm_b)
            do ind_i=1,nd_p1(me_b+1)
              i=n0_p1(me_b+1)+ind_i
              do l=0,lmpb(jsort)
                lp1=l+1
                do ie=1,ntle_pb(l,jsort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,jsort)
                    lm0=l+m+1
                    t0(lm0)=tmp(i,km)
                  enddo
                  call rotate1(t0,y0,l,u(1,ig),1)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,jsort)
                    lm0=l+m+1
                    f(ind_i,km,ind_r)=f(ind_i,km,ind_r)+y0(lm0)
                  enddo
                enddo
              enddo   !! over l
            enddo  !! over ind_i
          enddo    !! over ind_r
          deallocate(f0,tmp)
          isize=nd_p1(me_b+1)*n2pbm
          allocate(f0(isize,npnt,2))
          call dgemm('n','n',isize,npnt,ndim3_kk(me_k+1),1.d0,f,isize,
     &               trf(1,1,1),ndim3_kk(me_k+1),0.d0,f0(1,1,1),isize)
          call dgemm('n','n',isize,npnt,ndim3_kk(me_k+1),1.d0,f,isize,
     &               trf(1,1,2),ndim3_kk(me_k+1),0.d0,f0(1,1,2),isize)
          deallocate(f)
          if(nproc_k/=1) call dgop(f0,2*isize*npnt,'  +',comm_k)
          allocate(tc(n1pbm,n2pbm))
          iat=iat_sort(iatom,isort)
          jat=iat_sort(jatom,jsort)
          ind0=iopb(iat)-1
          jnd0=iopb(jat)-1
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            tc=(0.d0,0.d0)
            ij=0
            do j=1,n2pbm
              do ind_i=1,nd_p1(me_b+1)
                ij=ij+1
                i=n0_p1(me_b+1)+ind_i
                tc(i,j)=dcmplx(f0(ij,k,1),f0(ij,k,2))
              enddo
            enddo
            if(nproc_b/=1) call dgop(tc,2*n1pbm*n2pbm,'  +',comm_b)
            if(ncmpl==1) then
              km1=(jatom-1)*n2pbm
              do j=1,n2pbm
                km=(iatom-1)*n1pbm
                do i=1,n1pbm
                  p00(km+i,km1+j,ind_k)=p00(km+i,km1+j,ind_k)+tc(i,j)
                enddo
              enddo
            else
              call p_c_add_mm(pmm(1,ind_tau,ind_k),tc,n1pbm,n2pbm,ind0,
     &                        jnd0,n0,n1)
            endif
          enddo  !! over ind_k
          deallocate(tc,f0)
        enddo !! over iatom
      enddo   !! over jatom
      deallocate(nd_p1,n0_p1,trf)
      if(ncmpl==1) then
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          if(isort==jsort) then
            do j=1,nat2*n2pbm
              do i=1,j-1
                p00(j,i,ind_k)=conjg(p00(i,j,ind_k))
              enddo
            enddo
          endif
          call p_r_add_mm(pmm(1,ind_tau,ind_k),p00(1,1,ind_k),n1pbm,
     &                    n2pbm,n0,n1,nat1,nat2,isort,jsort,pnt(1,k))
        enddo
        deallocate(p00)
      endif
      end



      subroutine p_mm_sclp_fin(p_mm,nd_p2,n0_p2,n1pbm,n2pbm,nd_nrr,
     &                         n0_nrr,isort,jsort,nat1,nat2,nrr_reds,
     &                         r0_pairs,ig_pairs,ind_tau,pwsclp)
      use atom_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nd_p2,n0_p2,n1pbm,n2pbm,nd_nrr,n0_nrr,
     &                       isort,jsort,nat1,nat2,ind_tau,
     &                       nrr_reds(3,nqdiv*nat1*nat2),
     &                       r0_pairs(nqdiv,nat1,nat2),
     &                       ig_pairs(nqdiv,nat1,nat2)
      real*8, intent(in) :: p_mm(n1pbm,nd_p2,nd_nrr)
      real*8, intent(inout) :: pwsclp(ncmpl*nrow_pw_max,ncol_pw_max,
     &                                ndim3_tn,ndim3_k(me_k+1))
      integer :: j,ig,iatom,jatom,ir,ir0,ij,l,lp1,ie,m,lm,iat,ind0,
     &           km,lm0,i,nn0,nn1,ir1,ind_r,ind_i,ind_j,ind_k,isize,jat,
     &           k,jnd0,km1,n0,n1,na
      real*8 :: y0(2*maxpb+1),t0(2*maxpb+1),phase,pi2
      integer, allocatable :: nd_p1(:),n0_p1(:)
      real*8, allocatable :: tmp(:,:),t(:,:),f0(:,:,:),f(:,:,:),
     &                       trf(:,:,:)
      complex*16, allocatable :: tc(:,:),p00(:,:,:)
      pi2=pi+pi
      allocate(nd_p1(nproc_b))
      allocate(n0_p1(nproc_b))
      call size_shift_par(n1pbm,nproc_b,nd_p1,n0_p1)
      nn0=n0_nrr
      nn1=nn0+nd_nrr
      allocate(trf(ndim3_kk(me_k+1),npnt,2))
      do ind_r=1,ndim3_kk(me_k+1)
        ir=n3_mpi_kk(me_k+1)+ind_r
        do k=1,npnt
          phase=pi2*dot_product(pnt(:,k),rr_3(:,ir))
          trf(ind_r,k,1)=cos(phase)
          trf(ind_r,k,2)=-sin(phase)
        enddo
      enddo
c ------------------------------------------------------------------
      if(ncmpl==1) then
        allocate(p00(n1pbm*nat1,n2pbm*nat2,ndim3_k(me_k+1)))
        p00=(0.d0,0.d0)
      endif
      n0=n_mpi_pbmt(me_b+1)
      n1=n0+ndim_pbmt(me_b+1)
      do jatom=1,nat2
        na=nat1
        if(isort==jsort) na=jatom
        do iatom=1,na
          allocate(f0(n1pbm,nd_p2,ndim3_kk(me_k+1)))
          allocate(t(n1pbm,nqdiv))
          do j=1,nd_p2
            t=0.d0
            do ir=1,nqdiv
              ir0=r0_pairs(ir,iatom,jatom)
              if(ir0>nn0.and.ir0<=nn1) then
                ir1=ir0-nn0
                t(:,ir)=p_mm(:,j,ir1)
              endif
            enddo
            if(nproc_k/=1) call dgop(t,n1pbm*nqdiv,'  +',comm_k)
            do ind_r=1,ndim3_kk(me_k+1)
              ir=n3_mpi_kk(me_k+1)+ind_r
              f0(:,j,ind_r)=t(:,ir)
            enddo
          enddo
          deallocate(t)
          allocate(tmp(n1pbm,n2pbm))
          allocate(f(nd_p1(me_b+1),n2pbm,ndim3_kk(me_k+1)))
          f=0.d0
          do ind_r=1,ndim3_kk(me_k+1)
            ir=n3_mpi_kk(me_k+1)+ind_r
            ig=ig_pairs(ir,iatom,jatom)
            tmp=0.d0
            do ind_j=1,nd_p2
              j=n0_p2+ind_j
              do l=0,lmpb(isort)
                lp1=l+1
                do ie=1,ntle_pb(l,isort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,isort)
                    lm0=l+m+1
                    t0(lm0)=f0(km,ind_j,ind_r)
                  enddo
                  call rotate1(t0,y0,l,u(1,ig),1)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,isort)
                    lm0=l+m+1
                    tmp(km,j)=y0(lm0)
                  enddo
                enddo
              enddo   !! over l
            enddo   !! over ind_j
            if(nproc_b/=1) call dgop(tmp,n1pbm*n2pbm,'  +',comm_b)
            do ind_i=1,nd_p1(me_b+1)
              i=n0_p1(me_b+1)+ind_i
              do l=0,lmpb(jsort)
                lp1=l+1
                do ie=1,ntle_pb(l,jsort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,jsort)
                    lm0=l+m+1
                    t0(lm0)=tmp(i,km)
                  enddo
                  call rotate1(t0,y0,l,u(1,ig),1)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    km=ind_prod(ie,lm,jsort)
                    lm0=l+m+1
                    f(ind_i,km,ind_r)=f(ind_i,km,ind_r)+y0(lm0)
                  enddo
                enddo
              enddo   !! over l
            enddo  !! over ind_i
          enddo    !! over ind_r
          deallocate(f0,tmp)
          isize=nd_p1(me_b+1)*n2pbm
          allocate(f0(isize,npnt,2))
          call dgemm('n','n',isize,npnt,ndim3_kk(me_k+1),1.d0,f,isize,
     &               trf(1,1,1),ndim3_kk(me_k+1),0.d0,f0(1,1,1),isize)
          call dgemm('n','n',isize,npnt,ndim3_kk(me_k+1),1.d0,f,isize,
     &               trf(1,1,2),ndim3_kk(me_k+1),0.d0,f0(1,1,2),isize)
          deallocate(f)
          if(nproc_k/=1) call dgop(f0,2*isize*npnt,'  +',comm_k)
          allocate(tc(n1pbm,n2pbm))
          iat=iat_sort(iatom,isort)
          jat=iat_sort(jatom,jsort)
          ind0=iopb(iat)-1
          jnd0=iopb(jat)-1
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            tc=(0.d0,0.d0)
            ij=0
            do j=1,n2pbm
              do ind_i=1,nd_p1(me_b+1)
                ij=ij+1
                i=n0_p1(me_b+1)+ind_i
                tc(i,j)=dcmplx(f0(ij,k,1),f0(ij,k,2))
              enddo
            enddo
            if(nproc_b/=1) call dgop(tc,2*n1pbm*n2pbm,'  +',comm_b)
            if(ncmpl==1) then
              km1=(jatom-1)*n2pbm
              do j=1,n2pbm
                km=(iatom-1)*n1pbm
                do i=1,n1pbm
                  p00(km+i,km1+j,ind_k)=p00(km+i,km1+j,ind_k)+tc(i,j)
                enddo
              enddo
            else
              call p_c_add_mm_sclp(pwsclp(1,1,ind_tau,ind_k),tc,n1pbm,
     &                             n2pbm,ind0,jnd0,k,0)
              if(iat/=jat) then
                call p_c_add_mm_sclp(pwsclp(1,1,ind_tau,ind_k),tc,n1pbm,
     &                               n2pbm,ind0,jnd0,k,1)
              endif
            endif
          enddo  !! over ind_k
          deallocate(tc,f0)
        enddo !! over iatom
      enddo   !! over jatom
      deallocate(nd_p1,n0_p1,trf)
      if(ncmpl==1) then
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          if(isort==jsort) then
            do j=1,nat2*n2pbm
              do i=1,j-1
                p00(j,i,ind_k)=conjg(p00(i,j,ind_k))
              enddo
            enddo
          endif
          call p_r_add_mm_sclp(pwsclp(1,1,ind_tau,ind_k),p00(1,1,ind_k),
     &                         n1pbm,n2pbm,n0,n1,nat1,nat2,isort,jsort,
     &                         pnt(1,k),k)
        enddo
        deallocate(p00)
      endif
      end
