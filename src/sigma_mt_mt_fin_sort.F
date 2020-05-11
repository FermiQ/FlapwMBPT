      subroutine sigma_mt_mt_fin_sort(ind_k,k,ispin,isort,jsort,nfun1,
     &                                nfun2,pvt)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ind_k,k,ispin,isort,jsort,nfun1,nfun2
      complex*16, intent(inout) :: pvt(ndim_tb_nfun(me_t_b+1),nbndf)
      integer :: nz1,nz2,n,i,ind_i
      complex*16, allocatable :: tmp(:,:)
      n=n_bnd(k,ispin)
      nz1=io_lem(iat_1(isort))
      nz2=io_lem(iat_1(jsort))
      allocate(tmp(nfun1,n))
      call zgemm('n','n',nfun1,n,nfun2,(1.d0,0.d0),pvv(1,1,ind_k),nfun1,
     &           z_bnd(nz2,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp,nfun1)
      do ind_i=1,ndim_tb_nfun(me_t_b+1)
        i=nmpi_tb_nfun(me_t_b+1)+ind_i
        if(i<nz1.or.i>=nz1+nfun1) cycle
        pvt(ind_i,1:n)=pvt(ind_i,1:n)+tmp(i-nz1+1,:)
      enddo
      deallocate(tmp)
      end
