      subroutine g_rs_from_ks_sort(gr,gc,it,ind_tau,ispin,nf1,nf2,
     &                             nd_nrr,n0_nrr,nfun1,nfun2,isort,
     &                             jsort,nrr_maxs,nat1,nat2,nrr_reds)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ispin,nf1,nf2,nd_nrr,n0_nrr,nfun1,nfun2,
     &                       isort,jsort,nrr_maxs,nat1,nat2,it,
     &                       nrr_reds(3,nqdiv*nat1*nat2),ind_tau
      real*8, intent(out) :: gr(nf1,nf2,nd_nrr)
      complex*16, intent(out) :: gc(nf1,nf2,nd_nrr)
      integer :: k,k0,n,ind,j,ia,ib,ind_k,ig,iatom,ir,ir0,jatom,
     &           ind_ir0,ind_i,i,nz1,nz2,iat1,iat2,ind_r,i_tau
      real*8 :: pi2,tt(3),phase
      integer, allocatable :: nd_nf1(:),n0_nf1(:)
      real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp(:,:),g_tmp(:,:,:),tmp1(:,:),
     &                           tmp2(:,:),g_mm(:,:,:),t(:),
     &                           g1(:,:,:),tr1(:),tr2(:)
      allocate(g_tmp(nfun1,nfun2,ndim3_k(me_k+1)))
      pi2=pi+pi
      i_tau=me_t*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau
      nz1=io_lem(iat_1(isort))
      nz2=io_lem(iat_1(jsort))
      do ind_k=1,ndim3_k(me_k+1)
        k0=n3_mpi_k(me_k+1)+ind_k
        n=n_bnd(k0,ispin)
        allocate(tmp2(n,nfun2))
        if(ubi=='dft'.or.ubi==' hf') then
          allocate(gx_tau(n))
          call g_x_tau(ispin,k0,gx_tau,tau_mesh(i_tau),n,chem_pot)
          do j=1,nfun2
            do ib=1,n
              tmp2(ib,j)=gx_tau(ib)*conjg(z_bnd(nz2+j-1,ib,ind_k,ispin))
            enddo
          enddo
          deallocate(gx_tau)
        else
          allocate(tmp1(n,n))
          call ferm_unpack_tau(tmp1,g_full(1,1,1,ind_tau,ind_k,ispin),n,
     &                         n,nbndf,it)
          call zgemm('n','c',n,nfun2,n,(1.d0,0.d0),tmp1,n,
     &               z_bnd(nz2,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp2,n)
          deallocate(tmp1)
        endif
        call zgemm('n','n',nfun1,nfun2,n,(1.d0,0.d0),
     &             z_bnd(nz1,1,ind_k,ispin),nfun,tmp2,n,(0.d0,0.d0),
     &             g_tmp(1,1,ind_k),nfun1)
        deallocate(tmp2)
      enddo   !! over ind_k
      allocate(nd_nf1(nproc_b))
      allocate(n0_nf1(nproc_b))
      call size_shift_par(nfun1,nproc_b,nd_nf1,n0_nf1)
      allocate(g1(nd_nf1(me_b+1),nfun2,npnt))
      allocate(tmp1(nfun1,npnt))
      do j=1,nfun2
        tmp1=(0.d0,0.d0)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          tmp1(:,k)=g_tmp(:,j,ind_k)
        enddo
        if(nproc_k/=1) call dgop(tmp1,2*nfun1*npnt,'  +',comm_k)
        do ind_i=1,nd_nf1(me_b+1)
          i=n0_nf1(me_b+1)+ind_i
          g1(ind_i,j,:)=tmp1(i,:)
        enddo
      enddo
      deallocate(g_tmp,tmp1)
c ------------------------------------------------------------------
      allocate(tmp(nfun1,nfun2))
      allocate(tr1(nat1))
      allocate(tr2(nat2))
      allocate(g_mm(ndim3_kk(me_k+1),nd_nf1(me_b+1),nfun2))
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
        do ind_i=1,nd_nf1(me_b+1)
          i=n0_nf1(me_b+1)+ind_i
c ------------ Index from the right  ------------------------------
          call sym_val_left_sort(ig,g1(ind_i,:,k0),tmp(i,:),2,nfun2,
     &                           nat2,jsort,nf2)
          do iatom=1,nat2
            ia=nf2*(iatom-1)+1
            ib=ia+nf2-1
            tmp(i,ia:ib)=tmp(i,ia:ib)*conjg(tr2(iatom))
          enddo
        enddo
        if(nproc_b/=1) call dgop(tmp,2*nfun1*nfun2,'  +',comm_b)
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
          do ind_i=1,nd_nf1(me_b+1)
            i=n0_nf1(me_b+1)+ind_i
            g_mm(ind_k,ind_i,j)=t(i)
          enddo
        enddo
      enddo
      deallocate(tmp,g1,tr1,tr2,t)
c ------------------------------------------------------------------
      allocate(tmp(nqdiv,nd_nf1(me_b+1)))
      do j=1,nfun2
        tmp=(0.d0,0.d0)
        do ind_k=1,ndim3_kk(me_k+1)
          k=n3_mpi_kk(me_k+1)+ind_k
          call zone1_number(pnt(:,k),rb0,ndiv,ind)
          tmp(ind,:)=g_mm(ind_k,:,j)
        enddo
        if(nproc_k/=1) call dgop(tmp,2*nqdiv*nd_nf1(me_b+1),'  +',
     &                           comm_k)
        call fft3(ndiv(1),ndiv(2),ndiv(3),nd_nf1(me_b+1),tmp,1)
        do ind_r=1,ndim3_kk(me_k+1)
          ir=n3_mpi_kk(me_k+1)+ind_r
          g_mm(ind_r,:,j)=tmp(ir,:)
        enddo
      enddo
      deallocate(tmp)
      g_mm=g_mm/dfloat(nqdiv)
c ------------------------------------------------------------------
      allocate(tmp(nfun1,nqdiv))
      do j=1,nfun2
        tmp=(0.d0,0.d0)
        do ind_r=1,ndim3_kk(me_k+1)
          ir=n3_mpi_kk(me_k+1)+ind_r
          do ind_i=1,nd_nf1(me_b+1)
            i=n0_nf1(me_b+1)+ind_i
            tmp(i,ir)=g_mm(ind_r,ind_i,j)
          enddo
        enddo
        if(nproc_k/=1) call dgop(tmp,2*nqdiv*nfun1,'  +',comm_k)
        if(nproc_b/=1) call dgop(tmp,2*nqdiv*nfun1,'  +',comm_b)
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
      deallocate(g_mm,tmp,nd_nf1,n0_nf1)
      end
