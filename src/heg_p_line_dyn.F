      subroutine heg_p_line_dyn(spin_flips,dp)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      logical, intent(in) :: spin_flips
      real*8, intent(out) :: dp(n_line,0:n_nu)
      logical :: neg
      integer :: ind_nu,i_nu,i,nom,i_omega,nom_nu,it,iq,k,i_wi,ispin
      real*8 :: w1
      complex*16 :: cc
      real*8, allocatable :: w_n(:,:),
     &                       w_n_cheb(:,:),g_0(:,:,:,:),
     &                       g_0_cheb(:,:,:,:),pol(:,:),
     &                       conv_table(:,:,:)
      complex*16, allocatable :: g_om(:,:,:,:),g_om_x(:,:,:,:),
     &                           lambda_stat_nu_d(:,:,:,:),
     &                           lambda_stat_nu_d_sig2(:,:,:,:),
     &                           k_pw(:,:,:),lambda_dyn(:,:,:),
     &                           lambda_dyn_2(:,:,:)
	  allocate(lambda_stat_nu_d_sig2(nrr_div,nspin,ndim3_nu,n_line))
c ----------------------------------------------------------------
      allocate(conv_table(2:iter_ladder_xi,0:7,3))
	  allocate(lambda_stat_nu_d(nrr_div,nspin,ndim3_nu,n_line))
c --------------------------------------------------------------------
	  allocate(k_pw(nrr_div,0:n_tau,nspin))
c ------------------------------------------------------------------
      allocate(g_0(nstar,2,0:n_omega,nspin))
	  allocate(g_0_cheb(nstar,2,0:n_omega_asy-2,nspin))
	  do ispin=1,nspin
	    do i=1,nstar
	      do i_omega=0,n_omega
	        cc=g_k_omega_heg(i_omega,i,ispin)
     &	      -gx_k_omega_heg(i_omega,i,ispin)
            g_0(i,1,i_omega,ispin)=real(cc)
            g_0(i,2,i_omega,ispin)=imag(cc)
          enddo
        enddo
	    call dgemm('n','t',2*nstar,n_omega_asy-1,n_omega_asy-1,1.d0,
     &             g_0(1,1,n_omega_exa+1,ispin),2*nstar,cheb_omega(0,1),
     &             n_omega_asy-1,0.d0,g_0_cheb(1,1,0,ispin),2*nstar)
      enddo
c -----------------------------------------------------------------
      allocate(w_n(nstar,0:n_nu))
      allocate(w_n_cheb(nstar,0:n_nu_asy-2))
      if(spin_flips) then
        do i_nu=0,n_nu
          w_n(:,i_nu)=w_q_nu_heg(i_nu,:)
        enddo
	    call dgemm('n','t',nstar,n_nu_asy-1,n_nu_asy-1,1.d0,
     &             w_n(1,n_nu_exa+1),nstar,cheb_nu(0,1),n_nu_asy-1,
     &             0.d0,w_n_cheb,nstar)
      endif
c --------------------------------------------------------------------
      allocate(pol(n_line,0:n_nu))
      pol=0.d0
c --------------------------------------------------------------------
      if(iter_ladder_xi>1) conv_table=0.d0
	  do ind_nu=1,ndim3_nu
	    i_nu=me_t*ndim3_nu+ind_nu-1
	    nom=num_omega1_adapt(i_nu)
c --------------------------------------------------------------------
        w1=w_nu(i_nu)-pi/betta_t
        i_wi=1
        neg=.false.
        if(w1>0.d0) then
          neg=.true.
	      do i_omega=1,nom
            if(abs(w_om_adapt_nu(i_omega,i_nu)-w1)<1.d-7) i_wi=i_omega
          enddo
        endif
c --------------------------------------------------------------------
	    nom_nu=num_nu2_adapt(i_nu)
	    allocate(lambda_dyn(nrr_div,nom,nspin))
	    allocate(lambda_dyn_2(nrr_div,nom,nspin))
	    allocate(g_om(nstar,nom,2,nspin))
	    allocate(g_om_x(nstar,nom,2,nspin))
c ------ G-interpolation ------------------------------------------
        do ispin=1,nspin
	      do i_omega=1,nom
            w1=w_om_adapt_nu(i_omega,i_nu)
            call interp_g_omega(g_0(1,1,0,ispin),g_0_cheb(1,1,0,ispin),
     &                          g_om(1,i_omega,1,ispin),w1,0,nstar,1,1)
            do i=1,nstar
	          g_om_x(i,i_omega,1,ispin)=(1.d0,0.d0)
     & 	                       /dcmplx(chem_pot-e_star_x(i,ispin),w1)
            enddo
            w1=w_om_adapt_nu(i_omega,i_nu)-w_nu(i_nu)
            call interp_g_omega(g_0(1,1,0,ispin),g_0_cheb(1,1,0,ispin),
     &                          g_om(1,i_omega,2,ispin),w1,0,nstar,1,1)
            do i=1,nstar
	          g_om_x(i,i_omega,2,ispin)=(1.d0,0.d0)
     & 	                       /dcmplx(chem_pot-e_star_x(i,ispin),w1)
            enddo
          enddo
        enddo
        g_om=g_om+g_om_x
c ----------------------------------------------------------------
	    do iq=1,n_line
c -------------- Correction to the Polarizability ----------------------
          do ispin=1,nspin
            call k_heg(ispin,q_line(1,iq),nom,ind_nu,
     &                 lambda_dyn(1,1,ispin),k_pw(1,0,ispin),
     &                 lambda_stat_nu_d(1,ispin,ind_nu,iq),
     &                 g_om(1,1,1,ispin),g_om_x(1,1,1,ispin))
	        pol(iq,i_nu)=pol(iq,i_nu)-2.d0*real(k_pw(1,0,ispin))
	      enddo
        enddo   !! over iq
	    deallocate(lambda_dyn,g_om,g_om_x,lambda_dyn_2)
	  enddo   !! over ind_nu
	  if(nspin==1) pol=2.d0*pol
	  if(nproc_t/=1) then
	    call dgop(pol,n_line*(1+n_nu),'  +',comm_t)
        if(iter_ladder_xi>1) call dgop(conv_table,(iter_ladder_xi-1)*24,
     &                                 '  +',comm_t)
	  endif
	  if(maswrk.and.iter_ladder_xi>1) then
	    do k=1,3
	      write(iun,*)' Convergency of the ladder sequence, Q =',k
	      write(iun,'(a5,8(i3,6x))')'It\Nu',0,1,2,3,4,5,6,7
	      do it=2,iter_ladder_xi
	        write(iun,'(i3,8e9.2)')it,(conv_table(it,i_nu,k),i_nu=0,7)
          enddo
        enddo
      endif
c ------------------------------------------------------------------
      deallocate(k_pw,g_0,g_0_cheb,lambda_stat_nu_d,conv_table,w_n,
     &           w_n_cheb,lambda_stat_nu_d_sig2,pol)
      call timel('**** HEG_P_LINE_DYN finished *******')
      end
