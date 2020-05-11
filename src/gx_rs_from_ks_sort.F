      subroutine gx_rs_from_ks_sort(gr,gc,ispin,nf1,nf2,nd_nrr,n0_nrr,
     &                              nfun1,nfun2,isort,jsort,nrr_maxs,
     &                              nat1,nat2,nrr_reds)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ispin,nf1,nf2,nd_nrr,n0_nrr,nfun1,nfun2,
     &                       isort,jsort,nrr_maxs,nat1,nat2,
     &                       nrr_reds(3,nqdiv*nat1*nat2)
      real*8, intent(out) :: gr(nf1,nf2,nd_nrr)
      complex*16, intent(out) :: gc(nf1,nf2,nd_nrr)
      integer :: k,k0,nbnd,ind,j,ia,ib,ind_k,ig,iatom,ir,ir0,jatom,
     &           ind_ir0,ind_i,i,nz1,nz2,iat1,iat2,ind_r,npr
      real*8 :: pi2,tt(3),phase,gtm,de,green0
      integer, allocatable :: nd_nf1(:),n0_nf1(:)
      real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp(:,:),g_tmp(:,:,:),tmp1(:,:),
     &                           tmp2(:,:),g_mm(:,:,:),t(:),
     &                           g1(:,:,:),tr1(:),tr2(:)
      allocate(g_tmp(nfun1,nfun2,ndim3_k(me_k+1)))
      allocate(tmp2(nbndf,nfun2))
      pi2=pi+pi
      allocate(gx_tau(nbndf))
      g_tmp=(0.d0,0.d0)
      nz1=io_lem(iat_1(isort))
      nz2=io_lem(iat_1(jsort))
      do ind_k=1,ndim3_k(me_k+1)
        k0=n3_mpi_k(me_k+1)+ind_k
        nbnd=n_bnd(k0,ispin)
        if(ubi=='dft'.or.ubi==' hf') then
          do i=1,nbnd
            de=e_bnd(i,k0,ispin)-chem_pot
            gtm=-green0(de,betta_t)
            do j=1,nfun2
              tmp2(i,j)=gtm*conjg(z_bnd(nz2+j-1,i,ind_k,ispin))
            enddo
          enddo
        else
          call zgemm('n','c',nbnd,nfun2,nbnd,(1.d0,0.d0),
     &               g_full_0(1,1,ind_k,ispin),nbndf,
     &               z_bnd(nz2,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp2,
     &               nbndf)
        endif
        call zgemm('n','n',nfun1,nfun2,nbnd,(1.d0,0.d0),
     &             z_bnd(nz1,1,ind_k,ispin),nfun,tmp2,nbndf,(0.d0,0.d0),
     &             g_tmp(1,1,ind_k),nfun1)
      enddo   !! over ind_k
      deallocate(tmp2,gx_tau)
      npr=nproc_t*nproc_b
      allocate(nd_nf1(npr))
      allocate(n0_nf1(npr))
      call size_shift_par(nfun1,npr,nd_nf1,n0_nf1)
      if(nd_nf1(me_t_b+1)/=0) allocate(g1(nd_nf1(me_t_b+1),nfun2,npnt))
      allocate(tmp1(nfun1,npnt))
      do j=1,nfun2
        tmp1=(0.d0,0.d0)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          tmp1(:,k)=g_tmp(:,j,ind_k)
        enddo
        if(nproc_k/=1) call dgop(tmp1,2*nfun1*npnt,'  +',comm_k)
        if(nd_nf1(me_t_b+1)/=0) then
          do ind_i=1,nd_nf1(me_t_b+1)
            i=n0_nf1(me_t_b+1)+ind_i
            g1(ind_i,j,:)=tmp1(i,:)
          enddo
        endif
      enddo
      deallocate(g_tmp,tmp1)
c ------------------------------------------------------------------
      allocate(tmp(nfun1,nfun2))
      allocate(tr1(nat1))
      allocate(tr2(nat2))
      if(nd_nf1(me_t_b+1)/=0) allocate(g_mm(ndim3_kk(me_k+1),
     &                                      nd_nf1(me_t_b+1),nfun2))
      allocate(t(nfun1))
      do ind_k=1,ndim3_kk(me_k+1)
        k=n3_mpi_kk(me_k+1)+ind_k
        k0=i_kref(k)
        ig=k_group(k)
        do iat2=1,nat2
          iatom=iat_sort(iat2,jsort)
          tt=tshift(:,iatom,ig)
          phase=pi2*dot_product(pnt(:,k),tt)
          tr2(iat2)=dcmplx(cos(phase),sin(phase))
        enddo
        tmp=(0.d0,0.d0)
        if(nd_nf1(me_t_b+1)/=0) then
          do ind_i=1,nd_nf1(me_t_b+1)
            i=n0_nf1(me_t_b+1)+ind_i
c ------------ Index from the right  ------------------------------
            call sym_val_left_sort(ig,g1(ind_i,:,k0),tmp(i,:),2,nfun2,
     &                             nat2,jsort,nf2)
            do iatom=1,nat2
              ia=nf2*(iatom-1)+1
              ib=ia+nf2-1
              tmp(i,ia:ib)=tmp(i,ia:ib)*conjg(tr2(iatom))
            enddo
          enddo
        endif
        if(npr/=1) call dgop(tmp,2*nfun1*nfun2,'  +',comm_t_b)
        do iat1=1,nat1
          iatom=iat_sort(iat1,isort)
          tt=tshift(:,iatom,ig)
          phase=pi2*dot_product(pnt(:,k),tt)
          tr1(iat1)=dcmplx(cos(phase),sin(phase))
        enddo
        do j=1,nfun2
c ------------ Index from the left ------------------------------
          call sym_val_left_sort(ig,tmp(:,j),t,1,nfun1,nat1,isort,nf1)
          do iatom=1,nat1
            ia=nf1*(iatom-1)+1
            ib=ia+nf1-1
            t(ia:ib)=t(ia:ib)*tr1(iatom)
          enddo
          if(nd_nf1(me_t_b+1)/=0) then
            do ind_i=1,nd_nf1(me_t_b+1)
              i=n0_nf1(me_t_b+1)+ind_i
              g_mm(ind_k,ind_i,j)=t(i)
            enddo
          endif
        enddo
      enddo
      deallocate(tmp,tr1,tr2,t)
      if(nd_nf1(me_t_b+1)/=0) deallocate(g1)
c ------------------------------------------------------------------
      if(nd_nf1(me_t_b+1)/=0) then
        allocate(tmp(nqdiv,nd_nf1(me_t_b+1)))
        do j=1,nfun2
          tmp=(0.d0,0.d0)
          do ind_k=1,ndim3_kk(me_k+1)
            k=n3_mpi_kk(me_k+1)+ind_k
            call zone1_number(pnt(:,k),rb0,ndiv,ind)
            tmp(ind,:)=g_mm(ind_k,:,j)
          enddo
          if(nproc_k/=1) call dgop(tmp,2*nqdiv*nd_nf1(me_t_b+1),'  +',
     &                             comm_k)
          call fft3(ndiv(1),ndiv(2),ndiv(3),nd_nf1(me_t_b+1),tmp,1)
          do ind_r=1,ndim3_kk(me_k+1)
            ir=n3_mpi_kk(me_k+1)+ind_r
            g_mm(ind_r,:,j)=tmp(ir,:)
          enddo
        enddo
        deallocate(tmp)
        g_mm=g_mm/dfloat(nqdiv)
      endif
c ------------------------------------------------------------------
      allocate(tmp(nfun1,nqdiv))
      do j=1,nfun2
        tmp=(0.d0,0.d0)
        do ind_r=1,ndim3_kk(me_k+1)
          ir=n3_mpi_kk(me_k+1)+ind_r
          if(nd_nf1(me_t_b+1)/=0) then
            do ind_i=1,nd_nf1(me_t_b+1)
              i=n0_nf1(me_t_b+1)+ind_i
              tmp(i,ir)=g_mm(ind_r,ind_i,j)
            enddo
          endif
        enddo
        if(nproc_k/=1) call dgop(tmp,2*nqdiv*nfun1,'  +',comm_k)
        if(npr/=1) call dgop(tmp,2*nqdiv*nfun1,'  +',comm_t_b)
        do ind_ir0=1,nd_nrr
          ir0=n0_nrr+ind_ir0
          ir=nrr_reds(1,ir0)
          iatom=nrr_reds(2,ir0)
          jatom=nrr_reds(3,ir0)
          k=nf2*(jatom-1)
          if(j<=k.or.j>k+nf2) cycle
          k0=nf1*(iatom-1)
          if(irel<=1) then
            do i=k0+1,k0+nf1
              gr(i-k0,j-k,ind_ir0)=tmp(i,ir)
            enddo
          else if(irel==2) then
            do i=k0+1,k0+nf1
              gc(i-k0,j-k,ind_ir0)=tmp(i,ir)
            enddo
          endif
        enddo
      enddo
      if(nd_nf1(me_t_b+1)/=0) deallocate(g_mm)
      deallocate(tmp,nd_nf1,n0_nf1)
      end
