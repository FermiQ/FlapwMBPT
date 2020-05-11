      subroutine sigma_x_rs
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: ind_k,ispin,k
      real*8 :: dum
      complex*16 :: dumc
      real*8, allocatable :: g_r(:,:,:)
      complex*16, allocatable :: s_c1(:,:,:,:),g_c1(:,:,:,:)
      if(ubi/=' hf') call set_g_full_0
c -------- Precalculation of the matrix S*A ------------------------
      allocate(sa(nplw_gw,nbndf,ndim3_k(me_k+1),nrel,nspin))
      call sa_prep(sa)
c ------------------------------------------------------------------
c
c	---------- MT_MT Part ------------------------------------
c
      call sig_x_mm
      call timel('* Sigma_x_rs: MM finished **********')
c
c	---------- INT_MT Part ------------------------------------
c
      call v_im_r_from_k  !! here we open s_it_r0
      if(ndim_tb_nfun(me_t_b+1)/=0) then
        allocate(pvg(ndim_tb_nfun(me_t_b+1),nbasmpw,ndim3_k(me_k+1)))
      endif
      do ispin=1,nspin
        if(ndim_tb_nfun(me_t_b+1)/=0) pvg=(0.d0,0.d0)
        if(ndim_tb_nr_red(me_t_b+1)/=0) then
          if(irel<=1) then
            allocate(g_r(ndim_tb_nr_red(me_t_b+1),nfun,
     &                   ndim3_kk(me_k+1)))
            call gx_rs_mi_mem(g_r,dumc,ispin)
          else if(irel==2) then
            allocate(g_c1(ndim_tb_nr_red(me_t_b+1),nfun,
     &                    ndim3_kk(me_k+1),nrel))
            call gx_rs_mi_mem(dum,g_c1,ispin)
          endif
          allocate(s_c1(ndim_tb_nr_red(me_t_b+1),nfun,nrel,
     &                  ndim3_kk(me_k+1)))
          if(irel<=1) then
            call sigx_im_r(ispin,g_r,s_c1)
            deallocate(g_r)
          else if(irel==2) then
            call sigx_im_r_rel(ispin,g_c1,s_c1)
            deallocate(g_c1)
          endif
        else
          allocate(s_c1(1,1,1,1))
        endif
        call sigx_im_k_from_r(s_c1)
        deallocate(s_c1)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          call sigma_mt_int_fin(ind_k,k,ispin)
        enddo
      enddo
      if(ndim_tb_nfun(me_t_b+1)/=0) deallocate(pvg)
      if(ndim_tb_nr_red(me_t_b+1)/=0) deallocate(s_it_r0)
      call timel('* Sigma_x_rs: MI finished **********')
c
c	---------- INT_INT Part ------------------------------------
c
      allocate(pgg(nbasmpw,nbasmpw,ndim3_k(me_k+1)))
      do ispin=1,nspin
        pgg=(0.d0,0.d0)
        call sigx_ii_mem(ispin)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          call sigma_int_int_fin(ind_k,k,ispin)
        enddo
      enddo
      deallocate(pgg)
      deallocate(sa)
      if(ubi/=' hf') deallocate(g_full_0)
      call timel('* Sigma_x_rs: II finished **********')
      end
