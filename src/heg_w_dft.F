      subroutine heg_w_dft
	use atom_mod
	use etot_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
	implicit none
	integer :: i_nu,k
	real*8 :: d
      resp0_q_nu_heg=p_q_nu_heg
c ------- Exchange-Correlation Kernel ------------------------------
      f_xc=0.d0
      call heg_kernel_xc
c --------- Response_Enhanced ---------------------------------------
	do k=0,n_q_heg
	  do i_nu=0,n_nu
	    d=v_q_heg(k)+f_xc(k,i_nu)
	    resp_q_nu_heg(k,i_nu)=resp0_q_nu_heg(k,i_nu)
     &                         /(1.d0-resp0_q_nu_heg(k,i_nu)*d)
        enddo
      enddo
      call output_heg_f_xc
      call output_heg1_resp_nu(1)
      call heg_w_from_resp
      end