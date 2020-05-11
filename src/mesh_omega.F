      subroutine mesh_omega_spl
	use manager_mod
	use models_mod
	use parallel_mod
	use units_mod
	use vertex_mod
      implicit none
	integer :: i_omega,i,i0,j
	real*8 :: pi_bet,h1,fac,hx
	pi_bet=pi/betta_t
	h1=pi_bet+pi_bet
	allocate(w_omega(0:n_omega))
	allocate(mats_num_omega(0:n_omega))
	do i_omega=0,n_omega_exa
	  w_omega(i_omega)=pi_bet+i_omega*h1
	enddo
c ---------- Find asymptotic part of mesh ----------------------
      n_omega_maz=idnint((omega_max-pi_bet)/h1)
      omega_max=(2*n_omega_maz+1)*pi_bet
      w_omega(n_omega)=omega_max
      i=idnint((omega_geom-pi_bet)/h1)
      omega_geom=(2*i+1)*pi_bet
c --------- Geometric progression factor ------------------------
      hx=(omega_geom-w_omega(n_omega_exa))/h1/2.d0
      call geom_factor_finder(n_omega_geom,hx,fac)
      hx=h1*2
      do i=n_omega_exa+1,n_omega_exa+n_omega_geom
        w_omega(i)=w_omega(i-1)+hx
        j=idnint((w_omega(i)-pi_bet)/h1)
        w_omega(i)=(2*j+1)*pi_bet
        hx=hx*fac
      enddo
c -------------------------------------------------------------
      w_omega(n_omega-1)=omega_max-omega_last_step*h1
c -------------------------------------------------------------
      allocate(xm_omega(0:n_omega_asy+1))
      bb_omega=1.d0/w_omega(n_omega_exa+n_omega_geom)**2
      aa_omega=1.d0/omega_max**2-bb_omega
      xm_omega(0)=0.d0
      xm_omega(n_omega_asy+1)=1.d0
      xm_omega(n_omega_asy)=(1.d0/w_omega(n_omega-1)**2-bb_omega)
     &                     /aa_omega
      hx=xm_omega(n_omega_asy)/n_omega_asy
      do i=1,n_omega_asy-1
        i0=n_omega_exa+n_omega_geom+i
        xm_omega(i)=i*hx
        w_omega(i0)=1.d0/(aa_omega*xm_omega(i)+bb_omega)**(1.d0/2.d0)
        j=idnint((w_omega(i0)-pi_bet)/h1)
        w_omega(i0)=(2*j+1)*pi_bet
        xm_omega(i)=(1.d0/w_omega(i0)**2-bb_omega)/aa_omega
      enddo
c ------- Correspondence to true Matsubara's frequencies ---------- 
     	mats_num_omega(0)=0
	do i_omega=1,n_omega
	  mats_num_omega(i_omega)=nint((w_omega(i_omega)-pi_bet)/h1)
	enddo
	n_omega_maz=mats_num_omega(n_omega)
	if(maswrk)write(iun,*)
     &	 ' The number of Omega-points is ',n_omega+1
	if(maswrk)write(iun,*)
     &	 ' The number of Matsubara Omega-frequencies is ',n_omega_maz
c ---- Now - matrices for the transform OMEGA ----> TAU -------------
      allocate(tau_from_omega(0:n_omega,0:n_tau/2,2))
	call asympt_sums_omega
	call transf_tau_from_omega_spl
      allocate(dtau_from_omega(0:n_omega,2))
	call transf_dtau_from_omega_spl
      call timel('**** MESH_OMEGA_SPL finished *******')
      end