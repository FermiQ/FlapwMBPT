      subroutine sigx_mm_k_from_r_sort(sig_mm,nfun1,nfun2)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: nfun1,nfun2
      complex*16, intent(in) :: sig_mm(nfun1,nfun2,ndim3_kk(me_k+1))
      integer :: ind,ind_k,k,i,j,ind_ir,ir
      complex*16, allocatable :: p_c(:,:)
      allocate(p_c(nqdiv,nfun1))
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
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          call zone1_number(pnt(:,k),rb0,ndiv,ind)
          do i=1,nfun1
            pvv(i,j,ind_k)=p_c(ind,i)
          enddo
        enddo
      enddo
      deallocate(p_c)
      end
