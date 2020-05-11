      subroutine heg_susceptibility
      use atom_mod
      use etot_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer :: iq,iq0,i_nu
      real*8 :: d
      if(irel==2) return
c -------- LDA - based responses ----------------------------------
c -------- Save non-interacting response -------------------
      call heg1_pol
      resp0_q_nu_heg=p_q_nu_heg
      call output_heg_resp_nu(0)
      if(ubi=='dft'.and.chi_cmp/='000') then
c ------- Exchange-Correlation Kernel ------------------------------
        f_xc=0.d0
        f_xc_zz=0.d0
        call heg_kernel_xc
c --------- Response_Enhanced ---------------------------------------
	  do iq=1,n_line
	    iq0=kline_in_npnt(iq)
	    do i_nu=0,n_nu
	      d=(v_q_heg(iq0)+f_xc(iq0,i_nu))*resp0_q_nu_heg(i_nu,iq0)
	      resp_line_heg(iq,i_nu)=resp0_q_nu_heg(i_nu,iq0)/(1.d0-d)
	      d=f_xc_zz(iq0,i_nu)*resp0_q_nu_heg(i_nu,iq0)
	      resp_zz_line_heg(iq,i_nu)=resp0_q_nu_heg(i_nu,iq0)/(1.d0-d)
          enddo
        enddo
        call output_heg_resp_line
      endif
c ----------------------------------------------------------
      if(ubi/='dft') then
        if(chi_cmp(1:1)/='0'.or.chi_cmp(3:3)/='0') then
          if(vrt_x_appr(1:1)/='0') call response_uniform_heg
          if(vrt_x_appr(2:2)/='0') call response_static_heg
          if(vrt_x_appr(3:3)/='0') call response_dynamic_heg
        endif
      endif
      end