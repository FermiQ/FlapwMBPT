      subroutine p_r_tau_solid
	  use atom_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
      use scalapack_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer :: ispin
	  if(lsclp) then
        pw_sclp=0.d0
        do ispin=1,nspin
          call pol_sclp_spin_spin(ispin,ispin,pw_sclp)
        enddo  !! over ispin
c ---------------------------------------------------------------
        pw_sclp=pw_sclp*4.d0/nspin/nrel
	  else
	    p_w_q_tau_mm=0.d0
	    p_w_q_tau_mi=0.d0
	    p_w_q_tau_ii=0.d0
	    do ispin=1,nspin
	      call pol_0_spin_spin(ispin,ispin,p_w_q_tau_mm,p_w_q_tau_mi,
     &                         p_w_q_tau_ii)
        enddo  !! over ispin
c ---------------------------------------------------------------      
	    p_w_q_tau_mm=p_w_q_tau_mm*4.d0/nspin/nrel
	    p_w_q_tau_mi=p_w_q_tau_mi*4.d0/nspin/nrel
	    p_w_q_tau_ii=p_w_q_tau_ii*4.d0/nspin/nrel
	  endif
      call timel('** P_R_TAU_SOLID finished **********')
      end
