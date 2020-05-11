      subroutine self_energy
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8 :: conv_vh,conv_sx,conv_sc
      call set_adm(adm_sig)
      call self_energy_dft_basis
      if(ubi/=' hf') then
        call psi_sig_c_gw_g_full(psi_c_1_0)
        psi_c_1_0=psi_c_1_0/2.d0
      endif
      if(ubi=='psi') call add_sigma_vertex_bands(psi_sig)
      if(ubi=='bsp') call add_sigma_vertex_bands(psi2_sig)
c --------- Check of convergence -----------------------------------
      call conv_se(conv_vh,conv_sx,conv_sc)
      if(ubi==' hf') then
        if(maswrk) print '(a3,a9,e12.5,a8,e12.5)',
     &                    ubi,': CONV_VH',conv_vh,' CONV_SX',conv_sx
      else if(ubi/=' hf') then
        if(maswrk) print '(a3,a9,e12.5,a8,e12.5,a8,e12.5)',
     &                    ubi,': CONV_VH',conv_vh,' CONV_SX',conv_sx,
     &                                            ' CONV_SC',conv_sc
      endif
c --------------- Mixing the GW part of Sigma ----------------------
      call mixer_gw
c --------- Here we find a new spectrum ----------
      call bands_sig_0
      if(nbndf_bnd_0(2)/=0) call bands_vertex
      if(nbndf_bnd_0(2)/=0.and.mesh_k_c) then
c ---- Low energy bands for vertex corrections --------------
        call bands_vertex1
c ---- Integrals for PW-susceptibilities ---------------------------
        call psi_psi_pw_red
      endif
c ------------------------------------------------------------------
      v_h_bb_old=v_h_bb
      if(key_hybrid==1) v_xc_bb_old=v_xc_bb
      sigx_solid_k_old=sigx_solid_k
      if(ubi/=' hf') sig_c_tau_old=sig_c_tau
c ------------------------------------------------------------------
      call timel('**** SELF_ENERGY finished **********')
      end
