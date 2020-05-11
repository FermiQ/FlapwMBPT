      subroutine sigma_mt_int_fin(ind_k,k,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ind_k,k,ispin
      integer :: i,j,nd
      complex*16, allocatable :: sig_x1(:,:),tmp(:,:)
      allocate(sig_x1(nbndf,nbndf))
      sig_x1=(0.d0,0.d0)
      nd=ndim_tb_nfun(me_t_b+1)
      if(nd/=0) then
        allocate(tmp(nd,nbndf))
        call zgemm('n','n',nd,n_bnd(k,ispin),nbask(k),(1.d0,0.d0),
     &             pvg(1,1,ind_k),nd,ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &             (0.d0,0.d0),tmp,nd)
        call zgemm('c','n',n_bnd(k,ispin),n_bnd(k,ispin),nd,(1.d0,0.d0),
     &             z_bnd(nmpi_tb_nfun(me_t_b+1)+1,1,ind_k,ispin),nfun,
     &             tmp,nd,(0.d0,0.d0),sig_x1,nbndf)
        deallocate(tmp)
      endif
      if(nproc_t*nproc_b/=1) call dgop(sig_x1,2*nbndf**2,'  +',comm_t_b)
      do j=1,n_bnd(k,ispin)
        do i=1,j
          sig_x1(i,j)=sig_x1(i,j)+conjg(sig_x1(j,i))
        enddo
      enddo
      call pack_hrm_tb(sig_x1,sigx_solid_k(1,ind_k,ispin),
     &	               n_bnd(k,ispin),nbndf,1)
      deallocate(sig_x1)
      end
