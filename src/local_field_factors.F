	subroutine local_field_factors
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	implicit none
	integer :: i,i_nu
	real*8 :: xi0,p,v
	g_fac_s=0.d0
	g_fac_a=0.d0
	g_fac_n=0.d0
	f_xc=0.d0
	do i=2,nstar
	  v=v_q_heg(i)
	  do i_nu=0,n_nu
	    xi0=resp0_q_nu_heg(i_nu,i)
	    p=p_q_nu_heg(i_nu,i)
	    f_xc(i,i_nu)=(xi0-p)/p/xi0
	    g_fac_s(i,i_nu)=-f_xc(i,i_nu)/v
	  enddo
	enddo
	call output_heg_f_xc
	call output_heg_g_fac
      end
      
      
      subroutine local_field_factors_1
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	implicit none
	integer :: i,i_nu
	real*8 :: xi0,p,v
	g_fac_s=0.d0
	g_fac_a=0.d0
	g_fac_n=0.d0
	f_xc=0.d0
	do i=1,n_q_heg
	  v=v_q_heg(i)
	  do i_nu=0,n_nu
	    xi0=resp0_q_nu_heg(i,i_nu)
	    p=p_q_nu_heg(i,i_nu)
	    f_xc(i,i_nu)=(xi0-p)/p/xi0
	    g_fac_s(i,i_nu)=-f_xc(i,i_nu)/v
	  enddo
	enddo
	call output_heg1_f_xc
	call output_heg1_g_fac
	end