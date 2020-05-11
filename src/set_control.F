      subroutine set_control
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
      integer :: isort,l
      real*8 :: aa
      itest=0
      iprn=0
c ----------------------------------------------------------------
      n_omega=n_omega_exa+n_omega_asy+n_omega_geom+1
      n_nu=n_nu_exa+n_nu_asy+n_nu_geom+1
	  iter_max_scf=iter_dft+iter_hf+iter_gw+iter_qp+iter_psi+iter_bsp
	  if(iter_max_scf>iter_dft.and.rel_interst) then
	    if(maswrk) write(iun,*)' REL_INTERST=T for non DFT'
	    call ending
	  endif
      mesh_k_c=.false.
      if(vrt_x_appr/='000') mesh_k_c=.true.
      if(iter_psi+iter_bsp/=0) then
        if(psi_p/='00000'.or.bse_kernel_p/='0000') mesh_k_c=.true. 
        if(psi_sig/='00000'.or.psi2_sig/='00000') mesh_k_c=.true.
      endif
      if(iter_qp/=0.and.w_sc_qp=='non') then
        if(iter_gw/=0.and.w_sc_gw=='scf'.and.bse_kernel_p/='0000')
     &                   mesh_k_c=.true.
      endif
c ----------------------------------------------------------------
	  nqdiv=ndiv(1)*ndiv(2)*ndiv(3)
      if(b_extval>1.d-8) then
        if(irel<=1) then
          nspin_0=2
          nspin_1=2
          nspin=2
          magn=1
        endif
        if(irel==2) then
          nspin_0=1
          nspin_1=1
          nspin=1
          magn=2
        endif
      else
        aa=0.d0
        do isort=1,nsort
          aa=aa+abs(magn_shift(isort))
        enddo
        magn=1
        nspin_0=1
        if(aa>1.d-8) then
          if(irel<=1) nspin_0=2
          if(irel==2) magn=2
        endif
        nspin=1
        nspin_1=1
        if(irel<=1) then
          if(iter_h_ext/10000/=0) nspin=2
          if(iter_h_ext/10000/=0) nspin_1=2
        endif
      endif
c ------------------------------------------------------------------
      allocate(ri(0:40))
      allocate(rr(0:40))
      rr(0)=1.d0
      ri(0)=0.d0
      do l=1,40
         rr(l)=-ri(l-1)
         ri(l)=rr(l-1)
      enddo
c
c                 NLOC     ILDA
c
      nloc=iexch/100
      ilda=mod(iexch,100)
c
c                 NSPMAG
c
      nspmag=max(magn,nspin)
c -----------------------------------------------------------------
      betta_t=1.d0/temperature/boltzman
! ----------------------------------------------------------------
      omega_inf=1.d6
c ---------------------------------------------------------------      	
      l_double_freq=.false.	
      if(iter_max_scf/=iter_dft.and.vrt_x_appr/='000') 
     &                                      l_double_freq=.true.
      if(iter_psi+iter_bsp/=0.and.mesh_k_c) l_double_freq=.true.
      if(iter_hf/=0) then
        if(iter_gw/=0) then
          if(w_sc_gw=='non') l_double_freq=.true.
        else if(iter_qp/=0) then
          if(w_sc_qp=='non') l_double_freq=.true.
        else if(iter_psi/=0) then
          if(w_sc_psi=='non') l_double_freq=.true.
        else if(iter_bsp/=0) then
          if(w_sc_bsp=='non') l_double_freq=.true.
        endif
      endif
      if(iter_qp/=0.and.w_sc_qp=='non') then
        if(iter_gw/=0.and.w_sc_gw=='scf') l_double_freq=.true.
      endif
      if(l_double_freq) mesh_k_c=.true.
      z_fit=1
      emin_pdos_x=-5.d0
      emax_pdos_x=10.d0
      dp_tau_to_nu=0
      dsig_tau_to_omega=0
      l_grad_k=.false.
      if(maxval(n_k_int)/=1) l_grad_k=.true.
      if(hybrid_type==2) iexch=205
      chi_basis='PB'
      n_phead_pnt=2
      scale_psi=.false.
      end
