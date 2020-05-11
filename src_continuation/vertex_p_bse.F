      subroutine vertex_p_bse(key)
c     key = 0 - Based on DFT or HF Green's function
c     key = 1 - Based on QP Green's function
c     key = 2 - Based on general Green's function
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      integer :: ispin,mx2
      real*8, allocatable :: g_rr_tau(:,:,:,:,:,:),wst4loc(:,:,:,:,:)
      complex*16, allocatable :: tmatr(:,:,:,:,:,:)
      mx2=maxel_red**2
      if(bse_kernel_p=='2000') then
        call vertex_p_loc(iter_ladder_p,0)
      else if(bse_kernel_p=='2200') then
        call vertex_p_loc(iter_ladder_p,iter_ladder_p)
      else if(bse_kernel_p=='3000') then
        call vertex_p_vmt(iter_ladder_p,0,key)
        allocate(g_rr_tau(nfun_red,nfun_red,nqdiv_c,2,ndim3_tau,nspin))
        do ispin=1,nspin
          call g_real_space_mm(g_rr_tau(1,1,1,1,1,ispin),ispin)
        enddo
        allocate(wst4loc(mx2,mx2,natom,nspin,nspin))
        call w4_loc_stat(wst4loc)
        allocate(tmatr(n_lem2_red,n_lem2_red,ndim3_nu,ndimc_kk(me_k+1),
     &                 nspin,nspin))
        call t_matr_sf_gen(tmatr,g_rr_tau,wst4loc,'ph',iter_ladder_p+1)
        deallocate(g_rr_tau,wst4loc)
        call vertex_p_vmt_dir(tmatr)
        deallocate(tmatr)
      else if(bse_kernel_p=='3300') then
        call vertex_p_vmt(iter_ladder_p,iter_ladder_p,key)
      else if(bse_kernel_p=='4000') then
        call vertex_bnd_gw(iter_ladder_p,0,'wmt',key)
      else if(bse_kernel_p=='4400') then
        call vertex_bnd_gw(iter_ladder_p,iter_ladder_p,'wmt',key)
      else if(bse_kernel_p=='6000') then
        call vertex_bnd_gw(iter_ladder_p,0,'bnd',key)
      else if(bse_kernel_p=='6600') then
        call vertex_bnd_gw(iter_ladder_p,iter_ladder_p,'bnd',key)
      endif
      end
