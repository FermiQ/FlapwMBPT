      subroutine response_uniform_heg
      use atom_mod
      use etot_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8 :: d,dp
      resp_line_heg(1,0)=p_q_nu_heg(0,0)
      resp_zz_line_heg(1,0)=resp_line_heg(1,0)
c --- Vertex correction to the spin ZZ susceptibility --------------
c     if(ubi==' hf') call heg_p_line_uniform_hf(dp)
      if(ubi/=' hf') call heg_p_uniform(.false.,dp)
      resp_zz_line_heg(1,0)=resp_zz_line_heg(1,0)+dp !! Full ZZ response
c --- Vertex correction to the charge susceptibility ---------------
      if(ubi/=' hf') call heg_p_uniform(.true.,dp)
      dp=resp_line_heg(1,0)+dp  !! Full 00-polarizability
	d=dp*v_q_heg(1)
	resp_line_heg(1,0)=dp/(1.d0-d)
c ------------------------------------------------------------------
      if(maswrk) then
        write(iun,*)' Uniform 00-Polarizability:     ',dp
        write(iun,*)' Uniform spin ZZ-susceptibility:',
     &              resp_zz_line_heg(1,0)
      endif
      end