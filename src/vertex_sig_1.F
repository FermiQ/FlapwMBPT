      subroutine vertex_sig_1(psi,sig_1_1)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      character*1 :: psi
      complex*16, intent(out) :: sig_1_1(nbndf_bnd,nbndf_bnd,2,
     &                                   ndim3_tau,ndim_k_red(me_k+1),
     &                                   nspin)
      if(psi=='2') call vertex_sig_1_loc(sig_1_1)
      if(psi=='3') call vertex_sig_1_vmt(sig_1_1)
      if(psi=='4') call vertex_sig_1_wmt(sig_1_1)
      if(psi=='6') call vertex_sig_1_bnd(sig_1_1)
      call psi_sig_g_vrt(sig_1_1,psi_c_1_1)
      psi_c_1_1=psi_c_1_1/4.d0
      call psi_sig_vrt_g_full(sig_1_1,psi_c_1_1_fine)
      psi_c_1_1_fine=psi_c_1_1_fine/4.d0
      call timel('**** VERTEX_SIG_1 finished *********')
      end
