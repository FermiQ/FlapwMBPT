      subroutine tau_frequency_meshes_spl
c     Create TAU and frequency meshes -------------
      use manager_mod
      use models_mod
      use parallel_mod
      use units_mod
      use vertex_mod   
      implicit none
      real*8 :: pi_bet
      pi_bet=pi/betta_t
c ----------------- TAU mesh --------------------------------------- 
	  allocate(tau_mesh(0:n_tau))
	  call mesh_tau_spl
c ----------------- FREQUENCY meshes --------------------------------
      call mesh_nu_spl
c ------------------------------------------------------------------      
      call mesh_omega_spl
      if(l_double_freq) call omega_adapt_nu
c --- Prepare General transformations TAU ---> NU --------------------
      allocate(ab_nu_from_tau(0:n_tau/2+dp_tau_to_nu,0:n_nu,2))
      ab_nu_from_tau=0.d0
      call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),n_nu,w_nu(0),
     &                             ab_nu_from_tau,n_nu,dp_tau_to_nu)
c --- Basic AB transformations TAU ---> OMEGA -----------------------
      allocate(ab_omega_from_tau(0:n_tau/2+dsig_tau_to_omega,0:n_omega,
     &                           2))
      ab_omega_from_tau=0.d0
      call transf_ferm_from_tau_spl(n_tau/2,tau_mesh(0),n_omega,
     &                              w_omega(0),ab_omega_from_tau,
     &                              n_omega,dsig_tau_to_omega)
c ---------------------------------------------------------------
      if(l_double_freq) then
! -------- Here is adaptive w' mesh dense at 0 and w ----------------
!        call omega1_mesh_vertex
! -------- Here is adaptive v mesh dense at 0 and v ----------------
        call nu_adapt_nu
      endif
      if((iter_psi/=0.and.psi_sig/='00000').or.
     &   (iter_bsp/=0.and.psi2_sig/='00000')) then
        call omega_adapt_omega
        call nu_adapt_omega
      endif
      call timel('* TAU_FREQUNCY_MESHES_SPL finished *')
      end
