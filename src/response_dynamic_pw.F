      subroutine response_dynamic_pw(key)
!     key = 1 for 00 Polarizability
!     key = 2 for +- Polarizability
!     key = 3 for ZZ Polarizability
	  use atom_mod
	  use etot_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: key
      character*2 :: ch0
      integer :: iq,k,kc
      complex*16, allocatable :: pw(:,:),dpw(:),eps1(:),p0(:),eps1x(:)
      allocate(pw(0:n_nu,n_line))
      allocate(dpw(0:n_nu))
      allocate(p0(0:n_nu))
      allocate(eps1x(0:n_nu))
      allocate(eps1(0:n_nu))
      pw=(0.d0,0.d0)
      eps1=(0.d0,0.d0)
      if(key/=1) then
c --------- Non-Interacting case ------------------------------
        call resp_dynamic_prep_pw(pw,key)
        if(key==1) then
          call output_boson_q_nu(pw,'P0__00',n_line,'a')
        else if(key==2) then
          call output_boson_q_nu(pw,'P0__-+',n_line,'a')
        else if(key==3) then
          call output_boson_q_nu(pw,'P0__ZZ',n_line,'a')
        endif
c -------------- Real Axis --------------------------------------
        if(key==2) then
	      call boson_real_axis('P0__-+',pw,n_line,.false.,1)
          pw=conjg(pw)
	      call boson_real_axis('P0__+-',pw,n_line,.false.,1)
          pw=conjg(pw)
        else
          if(key==1) ch0='00'
          if(key==3) ch0='ZZ'
	      call boson_real_axis('P0__'//ch0,pw,n_line,.false.,2)
        endif
      endif
c --------- Vertex correction -------------------------------------
      dpw=(0.d0,0.d0)
      if(key==1) then
        call xi_eps1_00_pw_dynamic_gw(dpw,eps1,p0,eps1x)
      endif
      if(key==2) call vertex_xy_pw_dynamic_gw(dpw)
c     if(key==3) call vertex_zz_pw_dynamic_gw(dpw)
! ------------- Add non-interacting part ---------------
      do iq=1,n_line
	    k=kline_in_npnt(iq)
        kc=k_c_from_a(k)
        if(kc/=q_suscept) cycle
        dpw(:)=dpw(:)+pw(:,iq)
      enddo
! ----------------------------------------------------------      
      if(key==1) then
        if(p_reducible) then
          call output_boson_q_nu(dpw,'X___00',1,'c')
          call output_boson_q_nu(eps1,'E1__00',1,'c')
        else
          call output_boson_q_nu(dpw,'P___00',1,'c')
          call output_boson_q_nu(eps1,'E___00',1,'c')
        endif
        call output_boson_q_nu(p0,'P0C_00',1,'c')
        call output_boson_q_nu(eps1x,'E1X_00',1,'c')
      else if(key==2) then
        call output_boson_q_nu(dpw,'P1__+1',1,'c')
      else if(key==3) then
        call output_boson_q_nu(dpw,'P1__ZZ',1,'c')
      endif
c -------------- Real Axis --------------------------------------
      if(key==2) then
	    call boson_real_axis('P1__-+',dpw,1,.false.,1)
        dpw=conjg(dpw)
	    call boson_real_axis('P1__+-',dpw,1,.false.,1)
        dpw=conjg(dpw)
      else if(key==1) then
        if(p_reducible) then
	      call boson_real_axis('X_____',dpw,1,.false.,2)
	      call boson_real_axis('E1____',eps1,1,.true.,2)
        else
	      call boson_real_axis('P_____',dpw,1,.false.,2)
	      call boson_real_axis('E_____',eps1,1,.true.,2)
        endif
	    call boson_real_axis('E1X___',eps1x,1,.true.,2)
      else
        if(key==3) ch0='ZZ'
	    call boson_real_axis('P1__'//ch0,dpw,1,.true.,2)
      endif
      deallocate(pw,dpw,eps1,p0,eps1x)
      call timel('*** RESPONSE_DYNAMIC_PW finished ***')
      end
