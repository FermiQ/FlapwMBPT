      subroutine mesh_nu_spl
c ------- MESH for Nu ------------------------------------------------
	use manager_mod
	use models_mod
	use parallel_mod
	use units_mod
	use vertex_mod
      implicit none
	integer :: i_nu,i,i0,j,i_tau,i1_tau,ind_nu
	real*8 :: pi_bet,h1,fac,hx
	pi_bet=pi/betta_t
	h1=pi_bet+pi_bet
	allocate(w_nu(0:n_nu))
	allocate(mats_num_nu(0:n_nu))
	do i_nu=0,n_nu_exa
	  w_nu(i_nu)=i_nu*h1
	enddo
c ---------- Find asymptotic part of mesh ----------------------
      n_nu_maz=idnint(nu_max/h1)
      nu_max=n_nu_maz*h1
      w_nu(n_nu)=nu_max
      i=idnint(nu_geom/h1)
      nu_geom=i*h1
c --------- Geometric progression factor ------------------------
      hx=(nu_geom-w_nu(n_nu_exa))/h1/2.d0
      call geom_factor_finder(n_nu_geom,hx,fac)
      hx=h1*2
      do i=n_nu_exa+1,n_nu_exa+n_nu_geom
        w_nu(i)=w_nu(i-1)+hx
        j=idnint(w_nu(i)/h1)
        w_nu(i)=j*h1
        hx=hx*fac
      enddo
c -------------------------------------------------------------
      w_nu(n_nu-1)=nu_max-nu_last_step*h1
c -------------------------------------------------------------
      allocate(xm_nu(0:n_nu_asy+1))
      bb_nu=1.d0/w_nu(n_nu_exa+n_nu_geom)**2
      aa_nu=1.d0/nu_max**2-bb_nu
      xm_nu(0)=0.d0
      xm_nu(n_nu_asy+1)=1.d0
      xm_nu(n_nu_asy)=(1.d0/w_nu(n_nu-1)**2-bb_nu)/aa_nu
      hx=xm_nu(n_nu_asy)/n_nu_asy
      do i=1,n_nu_asy-1
        i0=n_nu_exa+n_nu_geom+i
        xm_nu(i)=i*hx
        w_nu(i0)=1.d0/sqrt(aa_nu*xm_nu(i)+bb_nu)
        j=idnint(w_nu(i0)/h1)
        w_nu(i0)=j*h1
        xm_nu(i)=(1.d0/w_nu(i0)**2-bb_nu)/aa_nu
      enddo
c ------- Correspondence to true Matsubara's frequencies ---------- 
     	mats_num_nu(0)=0
	do i_nu=1,n_nu
	  mats_num_nu(i_nu)=nint(w_nu(i_nu)/h1)
	enddo
	n_nu_maz=mats_num_nu(n_nu)
	if(maswrk)write(iun,*)
     &	 ' The number of Nu-points is ',n_nu+1
	if(maswrk)write(iun,*)
     &	 ' The number of Matsubara Nu-frequencies is ',n_nu_maz
c ---- Now - matrices for the transform NU ----> TAU -------------
      allocate(tau_from_nu(0:n_nu,0:n_tau/2,2))
	call asympt_sums_nu
	call transf_tau_from_nu_spl
c -------- For vertex correction to the Self Energy -------------
      allocate(tau_from_nu_d(ndim3_nu,0:n_tau,2))
	do i_tau=0,n_tau/2
	  i1_tau=n_tau-i_tau
	  do ind_nu=1,ndim3_nu
	    i_nu=me_t*ndim3_nu+ind_nu-1
	    tau_from_nu_d(ind_nu,i_tau,1)=tau_from_nu(i_nu,i_tau,1)/4
	    tau_from_nu_d(ind_nu,i_tau,2)=tau_from_nu(i_nu,i_tau,2)/4
	    tau_from_nu_d(ind_nu,i1_tau,1)=tau_from_nu(i_nu,i_tau,1)/4
	    tau_from_nu_d(ind_nu,i1_tau,2)=-tau_from_nu(i_nu,i_tau,2)/4
	  enddo
	enddo
c ---------- For Interpolation purposes -------------------------      
c      call for_nu_interpolation
      call timel('**** MESH_NU_SPL finished **********')
      end