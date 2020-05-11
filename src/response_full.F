      subroutine response_full
	  use atom_mod
	  use etot_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      complex*16, allocatable :: pw(:),e1w(:),rw(:),ew(:)
      allocate(pw(0:n_nu))
      allocate(rw(0:n_nu))
      allocate(ew(0:n_nu))
      allocate(e1w(0:n_nu))
      if(vrt_x_appr(2:2)=='1') call p00_1q_stat_w(pw,ew,e1w,rw)
      if(vrt_x_appr(3:3)=='1') call p00_1q_dyn_w(pw,ew,e1w,rw)
      call output_boson_q_nu(rw,'X___00',1,'c')
      call output_boson_q_nu(e1w,'E1__00',1,'c')
      call output_boson_q_nu(pw,'P___00',1,'c')
      call output_boson_q_nu(ew,'E___00',1,'c')
c -------------- Real Axis --------------------------------------
      call boson_real_axis('X_____',rw,1,.false.,2)
	  call boson_real_axis('E1____',e1w,1,.true.,2)
      call boson_real_axis('P_____',pw,1,.false.,2)
	  call boson_real_axis('E_____',ew,1,.true.,2)
      deallocate(pw,rw,e1w,ew)
      call timel('*** RESPONSE_DYNAMIC_FULL finished *')
      end
