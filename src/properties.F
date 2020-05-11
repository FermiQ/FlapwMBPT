      subroutine properties(iexit1)
      use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
      implicit none
	  integer, intent(inout) :: iexit1
	  if(iter<=iter_dft+iter_hf.or.restart_begin<=restart_end) then
	    call output_band_energy
c ------------ Numbers of Occupation -------------------------------
	    call occ_number
c ------ Output GW Energy ------------------------------------------
	    if(ubi/=' qp') call etot_gw_2
      endif
	  if(iexit1==1) then
	    if(iter<=iter_dft+iter_hf.or.restart_begin<=restart_end) then
c          if(ubi/='dft'.and.ubi/=' hf') call local_theta_kernel
c -------------------- Bands plotting ------------------------------
          if(ubi=='dft') call plot_bands_lapw !! DFT bands
          if(ubi=='dft'.or.ubi==' hf') then
            call bands_interp(e_bnd)
            call output_boltztrap(e_bnd)
          endif
          if(ubi/='dft'.and.ubi/=' hf') then
            call bands_interp(e_qp)
            call output_boltztrap(e_qp)
          endif
c ------------------------------------------------------------------
c          call forces
        endif
c        call optics
        if(ubi=='dft') then
          call susceptibility_dft
        else if(ubi==' hf') then
          call susceptibility_hf
        else
          call susceptibility_gw
        endif
c        if(ubi/='dft'.and.ubi/=' hf') call local_theta_kernel
	  endif
      if(ubi==' qp'.or.ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
        if(iter==restart_end.and.restart_begin>restart_end) then
          sig_c_tau=sig_c_tau_old
        endif
   	    call output_dos
        call output_pdos_x
      else if(iexit1==1) then
        call output_dos
        call output_pdos_x
      endif
      call timel('**** PROPERTIES finished ***********')
      end
