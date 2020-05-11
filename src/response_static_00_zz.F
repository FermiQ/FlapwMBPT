      subroutine response_static_00_zz
	  use atom_mod
	  use etot_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer :: ind_k,k,iq,n_pbt,kc,ka
      real*8, allocatable :: pw_00(:),pw_zz(:),dpw_00(:),dpw_zz(:),
     &                       v_tmpr(:,:)
      complex*16, allocatable :: p_non(:,:,:,:),pw(:,:),
     &                           vc(:,:,:),v_tmp(:,:),dpb_00(:,:,:),
     &                           dpb_0z(:,:,:),dpb_z0(:,:,:),
     &                           dpb_zz(:,:,:),tmp(:,:)
      allocate(p_non(n_pbtot,n_pbtot,nc_line,3))
      p_non=(0.d0,0.d0)
      call resp_static_prep_00_zz(p_non)
      allocate(pw_00(nc_line))
      allocate(pw_zz(nc_line))
      allocate(dpw_00(nc_line))
      allocate(dpw_zz(nc_line))
      allocate(pw(n_pbtot,nc_line))
      allocate(vc(n_pbtot,n_pbtot,nc_line))
      pw=(0.d0,0.d0)
      vc=(0.d0,0.d0)
      pw_00=0.d0
      pw_zz=0.d0
      dpw_00=0.d0
      dpw_zz=0.d0
	  do ind_k=1,ndim3_k(me_k+1)
	    k=n3_mpi_k(me_k+1)+ind_k
        kc=k_c_from_a(k)
        if(kc==0) cycle
	    iq=k_npnt_in_line_c(kc)
        if(iq==0) cycle
	    n_pbt=n_pbmt+nplwgw(k)
	    allocate(v_tmp(n_pbt,n_pbt))
	    pw(1:n_pbt,iq)=pw_pb(1:n_pbt,ind_k)
        if(ncmpl==1) then
          allocate(v_tmpr(n_pbt,n_pbt))
          call v_coul_full_r(nplwgw(k),v_tmpr,ind_k,0)
          call pb_c_from_r(n_pbt,n_pbt,n_pbt,v_tmp,v_tmpr,pnt(1,k))
          deallocate(v_tmpr)
        else
          call v_coul_full(nplwgw(k),v_tmp,ind_k,0)
        endif
	    vc(1:n_pbt,1:n_pbt,iq)=v_tmp
	    deallocate(v_tmp)
      enddo
	  if(nproc_k/=1) then
	    call dgop(pw,2*n_pbtot*nc_line,'  +',comm_k)
	    call dgop(vc,2*n_pbtot**2*nc_line,'  +',comm_k)
	  endif
c --------- Non-Interacting case ------------------------------
	  do iq=1,nc_line
	    k=kline_in_npnt_c(iq)
        ka=k_a_from_c(k)
	    n_pbt=n_pbmt+nplwgw(ka)
        call pw_from_pb_line(n_pbt,p_non(1,1,iq,1),pw_00(iq),3,pw(1,iq))
      enddo
      if(chi_cmp(3:3)/='0') pw_zz=pw_00
      if(chi_cmp(1:1)=='0') pw_00=0.d0
      if(nspin==1) then
        call output_resp_mbpt_static(p_non,p_non,p_non,pw_00,pw_zz,vc,
     &                               pw,0)
      else if(nspin==2) then
        call output_resp_mbpt_static(p_non(1,1,1,1),p_non(1,1,1,2),
     &                               p_non(1,1,1,2),pw_00,pw_zz,vc,pw,
     &                               0)
      endif
c --------- Leading term only vertex correction ---------------
      if(chi_cmp(1:1)/='0') call vertex_pol_pw_static_gw(0,dpw_00)
      if(chi_cmp(3:3)/='0') call vertex_pol_pw_static_gw(1,dpw_zz)
      if(chi_cmp(1:1)/='0') pw_00=pw_00+dpw_00
      if(chi_cmp(3:3)/='0') pw_zz=pw_zz+dpw_zz
c --------- Full vertex correction ----------------------------
	  allocate(dpb_00(n_pbtot_red,n_pbtot_red,nc_line))
	  allocate(dpb_0z(n_pbtot_red,n_pbtot_red,nc_line))
	  allocate(dpb_z0(n_pbtot_red,n_pbtot_red,nc_line))
	  allocate(dpb_zz(n_pbtot_red,n_pbtot_red,nc_line))
      dpb_00=(0.d0,0.d0)
      dpb_zz=(0.d0,0.d0)
      allocate(tmp(n_pbtot,n_pbtot))
      if(chi_cmp(1:1)/='0') 
     &    call vertex_pol_pb_static_gw(0,dpb_00,dpb_z0)   ! 00-vertex
      if(nspin==2.and.chi_cmp(3:3)/='0') 
     &    call vertex_pol_pb_static_gw(1,dpb_zz,dpb_0z)   ! ZZ-vertex
      if(chi_cmp(1:1)/='0') then
	    do iq=1,nc_line
	      k=kline_in_npnt_c(iq)
          ka=k_a_from_c(k)
	      n_pbt=n_pbmt+nplwgw(ka)
	      call fpb_from_rpb_full(dpb_00(1,1,iq),tmp,k,n_pbtot_red,
     &	                         n_pbtot,0)
	      p_non(1:n_pbt,1:n_pbt,iq,1)=p_non(1:n_pbt,1:n_pbt,iq,1)
     &                               +tmp(1:n_pbt,1:n_pbt)
        enddo
      endif
      if(nspin==2.and.chi_cmp(3:3)/='0') then
	    do iq=1,nc_line
	      k=kline_in_npnt_c(iq)
          ka=k_a_from_c(k)
	      n_pbt=n_pbmt+nplwgw(ka)
	      call fpb_from_rpb_full(dpb_z0(1,1,iq),tmp,k,n_pbtot_red,
     &	                         n_pbtot,0)
	      p_non(1:n_pbt,1:n_pbt,iq,3)=p_non(1:n_pbt,1:n_pbt,iq,2)
     &                               +tmp(1:n_pbt,1:n_pbt)
	      call fpb_from_rpb_full(dpb_0z(1,1,iq),tmp,k,n_pbtot_red,
     &	                         n_pbtot,0)
	      p_non(1:n_pbt,1:n_pbt,iq,2)=p_non(1:n_pbt,1:n_pbt,iq,2)
     &                               +tmp(1:n_pbt,1:n_pbt)
        enddo
      endif
      call output_p_static_q(dpb_00,dpb_zz,dpw_00,dpw_zz)
      call output_resp_mbpt_static(p_non(1,1,1,1),p_non(1,1,1,2),
     &                             p_non(1,1,1,3),pw_00,pw_zz,vc,pw,1)
      deallocate(p_non,pw_00,pw_zz,pw,vc,dpw_00,dpw_zz,dpb_00,
     &           dpb_0z,dpb_z0,dpb_zz,tmp)
      call timel('*** RESPONSE_STATIC_00_ZZ finished *')
      end
