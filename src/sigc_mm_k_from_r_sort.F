      subroutine sigc_mm_k_from_r_sort(ispin,it,sig_mm,nfun1,
     &                                 nfun2,pvt,isort,jsort)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ispin,it,nfun1,nfun2,isort,jsort
      complex*16, intent(in) :: sig_mm(nfun1,nfun2,ndim3_kk(me_k+1))
      complex*16, intent(inout) :: pvt(ndim_b_nfun(me_b+1),nbndf,2,
     &                                 ndim3_k(me_k+1))
      integer :: ind,ind_k,k,n,j,ind_ir,ir,i,nz1,nz2,ind_i
      complex*16, allocatable :: p_c(:,:),pcc(:,:),pc1(:,:,:),pc2(:,:)
      allocate(p_c(nqdiv,nfun1))
      allocate(pcc(nfun1,nqdiv))
      allocate(pc1(nfun1,nfun2,ndim3_k(me_k+1)))
      pc1=(0.d0,0.d0)
      do j=1,nfun2
        p_c=(0.d0,0.d0)
        do ind_ir=1,ndim3_kk(me_k+1)
          ir=n3_mpi_kk(me_k+1)+ind_ir
          do i=1,nfun1
            p_c(ir,i)=sig_mm(i,j,ind_ir)
          enddo
        enddo
        if(nproc_k/=1) call dgop(p_c,2*nqdiv*nfun1,'  +',comm_k)
        call fft3(ndiv(1),ndiv(2),ndiv(3),nfun1,p_c,-1)
        do k=1,nqdiv
          call zone1_number(pnt(:,k),rb0,ndiv,ind)
          do i=1,nfun1
            pcc(i,k)=p_c(ind,i)
          enddo
        enddo  !! over k
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          pc1(:,j,ind_k)=pc1(:,j,ind_k)+pcc(:,k)
        enddo  !! over ind_k
      enddo  !! over j
      deallocate(p_c,pcc)
c ------------------------------------------------------------------
      nz1=io_lem(iat_1(isort))
      nz2=io_lem(iat_1(jsort))
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        n=n_bnd(k,ispin)
        allocate(pc2(nfun1,n))
        call zgemm('n','n',nfun1,n,nfun2,(1.d0,0.d0),pc1(1,1,ind_k),
     &             nfun1,z_bnd(nz2,1,ind_k,ispin),nfun,(0.d0,0.d0),pc2,
     &             nfun1)
        do ind_i=1,ndim_b_nfun(me_b+1)
          i=nmpi_b_nfun(me_b+1)+ind_i
          if(i<nz1.or.i>=nz1+nfun1) cycle
          pvt(ind_i,1:n,it,ind_k)=pvt(ind_i,1:n,it,ind_k)+pc2(i-nz1+1,:)
        enddo
        deallocate(pc2)
      enddo  !! over ind_k
      deallocate(pc1)
      end
