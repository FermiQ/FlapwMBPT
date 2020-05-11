      subroutine ws_k_nu_solid
	  use atom_mod
	  use etot_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  logical :: w_last,w_red
      integer :: k,ind_k,npbt1,it0
      real*8, allocatable :: w_out_dir(:,:),e_out_dir(:,:),
     &                       e1_out_dir(:,:)
      complex*16, allocatable :: eout(:,:,:),e1out(:,:,:)
	  w_last=.false.
	  it0=iter_dft+iter_hf+iter_gw
	  if(ubi==' gw'.and.iter==it0) w_last=.true.
	  if(ubi==' qp'.and.iter==it0+iter_qp) w_last=.true.
	  if(ubi==' gw'.and.w_sc_gw=='non') w_last=.true.
	  if(ubi==' qp'.and.w_sc_qp=='non') w_last=.true.
	  if(ubi=='psi'.and.iter==it0+iter_qp+iter_psi) w_last=.true.
	  if(ubi=='bsp'.and.iter==it0+iter_qp+iter_psi+iter_bsp)
     &                                              w_last=.true.
	  if(ubi=='psi'.and.w_sc_psi=='non') w_last=.true.
	  if(ubi=='bsp'.and.w_sc_bsp=='non') w_last=.true.
      w_red=.false.
      if(nbndf_bnd_0(2)/=0) then
        if(ubi==' gw'.and.w_last) then
          if(iter_psi+iter_bsp/=0) w_red=.true.
          if(iter_ladder_xi/=0) w_red=.true.
          if(iter_qp/=0.and.w_sc_qp=='non') w_red=.true.
        endif
        if(ubi==' qp'.and.w_last) then
          if(iter_psi+iter_bsp/=0) w_red=.true.
          if(iter_ladder_xi/=0) w_red=.true.
        endif
        if(ubi=='psi') then
          if(w_sc_psi=='scf') w_red=.true.
          if(w_sc_psi/='scf'.and.psi_sig/='00000') w_red=.true.
        endif
        if(ubi=='bsp') then
          if(w_sc_bsp=='scf') w_red=.true.
          if(w_sc_bsp/='scf'.and.psi_sig/='00000') w_red=.true.
        endif
      endif
      w_red_print=w_red
	  allocate(e_output(0:n_nu,n_line))
	  e_output=0.d0
	  allocate(w_output(0:n_nu,n_line))
	  w_output=0.d0
	  allocate(p_output(0:n_nu,n_line))
	  p_output=0.d0
	  if(mode_bts==1) then
        allocate(e_out_dir(n_phead_dir,ndim3_nu))
        e_out_dir=0.d0
        allocate(e1_out_dir(n_phead_dir,ndim3_nu))
        e1_out_dir=0.d0
        allocate(w_out_dir(n_phead_dir,ndim3_nu))
        w_out_dir=0.d0
	  endif
	  npbt1=n_pbmt+nplwgw(1)
c ------------------------------------------------------------------
	  e_c_solid_pw=0.d0
	  ws_head=0.d0
	  ws_head_nu=0.d0
	  phi_functional_c=0.d0
      if(w_red) then
        w_nu_q_b=0.d0
        if(.not.allocated(w11_red))
     &    allocate(w11_red(n_opt_pb_red,n_opt_pb_red,ndim3_nu))
        w11_red=(0.d0,0.d0)
      endif
      if(ncmpl==1) then
        if(.not.allocated(ws_head_tens_r))
     &                     allocate(ws_head_tens_r(3,3,ndim3_nu))
        ws_head_tens_r=0.d0
        if(.not.allocated(ws_wing_r))
     &                     allocate(ws_wing_r(3,n_opt_pb,ndim3_nu))
        ws_wing_r=0.d0
        if(.not.allocated(ws01_met_r)) allocate(ws01_met_r(n_opt_pb))
        ws01_met_r=0.d0
      else
        if(.not.allocated(ws_head_tens))
     &                     allocate(ws_head_tens(3,3,ndim3_nu))
        ws_head_tens=(0.d0,0.d0)
        if(.not.allocated(ws_wing))
     &                     allocate(ws_wing(3,n_opt_pb,ndim3_nu))
        ws_wing=(0.d0,0.d0)
        if(.not.allocated(ws01_met)) allocate(ws01_met(n_opt_pb))
        ws01_met=(0.d0,0.d0)
      endif
      allocate(eout(3,3,0:n_nu))
      allocate(e1out(3,3,0:n_nu))
      eout=(0.d0,0.d0)
      e1out=(0.d0,0.d0)
	  do ind_k=1,ndim3_k(me_k+1)
	    k=n3_mpi_k(me_k+1)+ind_k
	    if(k==1) then
	      if(mode_bts==0) then
            call ws0_q0(w_last,w_red,eout,e1out)
          else
            call ws0_q0_fxc(w_out_dir,e_out_dir,e1_out_dir)
          endif
        else
          if(mode_bts==0) then
            call ws_k_nu_solid_0(ind_k,w_last,w_red)
          else
            call ws_q_fxc(ind_k)
          endif
	    endif
      enddo  !! over ind_k
c ------------------------------------------------------------------
	  if(nproc_k/=1) then
	    call dgop(e_c_solid_pw,1,'  +',comm_k)
	    call dgop(phi_functional_c,1,'  +',comm_k)
	    call dgop(w_output,n_line*(n_nu+1),'  +',comm_k)
	    call dgop(p_output,n_line*(n_nu+1),'  +',comm_k)
        call dgop(ws_head,n_tau/2+1,'  +',comm_k)
        call dgop(ws_head_nu,n_nu+1,'  +',comm_k)
        if(ncmpl==1) then
          call dgop(ws_head_tens_r,9*ndim3_nu,'  +',comm_k)
          call dgop(ws_wing_r,3*n_opt_pb*ndim3_nu,'  +',comm_k)
          call dgop(ws01_met_r,n_opt_pb,'  +',comm_k)
        else
          call dgop(ws_head_tens,18*ndim3_nu,'  +',comm_k)
          call dgop(ws_wing,6*n_opt_pb*ndim3_nu,'  +',comm_k)
          call dgop(ws01_met,2*n_opt_pb,'  +',comm_k)
        endif
        call dgop(eout,18*(n_nu+1),'  +',comm_k)
        call dgop(e1out,18*(n_nu+1),'  +',comm_k)
	    call dgop(e_output,n_line*(n_nu+1),'  +',comm_k)
        if(w_red) call dgop(w11_red,2*n_opt_pb_red**2*ndim3_nu,'  +',
     &                      comm_k)
	  endif
	  if(nproc_t/=1) then
	    call dgop(e_c_solid_pw,1,'  +',comm_t)
	    call dgop(phi_functional_c,1,'  +',comm_t)
	    call dgop(w_output,n_line*(n_nu+1),'  +',comm_t)
	    call dgop(p_output,n_line*(n_nu+1),'  +',comm_t)
        call dgop(eout,18*(n_nu+1),'  +',comm_t)
        call dgop(e1out,18*(n_nu+1),'  +',comm_t)
	    call dgop(e_output,n_line*(n_nu+1),'  +',comm_t)
	  endif
c -------------------------------------------------------------- 
      call output_e_tensor(eout,'_E')
      call output_e_tensor(e1out,'E1')
      deallocate(eout,e1out)
	  call output_w
	  deallocate(w_output,p_output,e_output)
	  if(mode_bts==1) then
	    call output_e_dir(ph_dir,w_out_dir,e_out_dir,e1_out_dir)
	    deallocate(e_out_dir,e1_out_dir,w_out_dir)
	  endif
      end
