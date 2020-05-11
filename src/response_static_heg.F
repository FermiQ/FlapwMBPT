      subroutine response_static_heg
      use atom_mod
      use etot_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer :: iq,iq0
      real*8 :: d
      real*8, allocatable :: dp(:)
	  do iq=1,n_line
	    iq0=kline_in_npnt(iq)
        resp_line_heg(iq,0)=p_q_nu_heg(0,iq0)
      enddo
      resp_zz_line_heg(:,0)=resp_line_heg(:,0)
c ------------------------------------------------------------------
      allocate(dp(n_line))
c --- Vertex correction to the spin ZZ susceptibility --------------
      if(ubi/=' hf') call heg_p_line_stat(.false.,dp)
      resp_zz_line_heg(:,0)=resp_zz_line_heg(:,0)+dp !! Full ZZ response
c --- Vertex correction to the charge susceptibility ---------------
      if(ubi/=' hf') call heg_p_line_stat(.true.,dp)
      dp=resp_line_heg(:,0)+dp  !! Full 00-polarizability
	  do iq=1,n_line
	    iq0=kline_in_npnt(iq)
	    d=dp(iq)*v_q_heg(iq0)
	    resp_line_heg(iq,0)=dp(iq)/(1.d0-d)
      enddo
      deallocate(dp)
c ------------------------------------------------------------------
      call output_heg_resp_static
      end
