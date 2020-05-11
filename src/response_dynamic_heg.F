      subroutine response_dynamic_heg
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
      real*8, allocatable :: dp(:,:)
	do iq=1,n_line
	  iq0=kline_in_npnt(iq)
        resp_line_heg(iq,:)=p_q_nu_heg(:,iq0)
      enddo
      resp_zz_line_heg=resp_line_heg
c ------------------------------------------------------------------
      allocate(dp(n_line,0:n_nu))
c --- Vertex correction to the spin ZZ susceptibility --------------
      if(ubi/=' hf') call heg_p_line_dyn(.false.,dp)
      resp_zz_line_heg=resp_zz_line_heg+dp  !! Full ZZ response
c --- Vertex correction to the charge susceptibility ---------------
      if(ubi/=' hf') call heg_p_line_dyn(.true.,dp)
      dp=resp_line_heg+dp  !! Full 00-polarizability
	do iq=1,n_line
	  iq0=kline_in_npnt(iq)
	  do i_nu=0,n_nu
	    d=dp(iq,i_nu)*v_q_heg(iq0)
	    resp_line_heg(iq,i_nu)=dp(iq,i_nu)/(1.d0-d)
        enddo
      enddo
      deallocate(dp)
c ------------------------------------------------------------------
      call output_heg_resp_line
      end