      subroutine transf_tau_from_omega_spl
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: i_omega,i_tau
      real*8 :: c1,w
      real*8, allocatable :: tmp(:,:,:)
c ----------------------------------------------------------------
      tau_from_omega=0.d0
c ---------- Exact contribution ----------------------------------  
      do i_omega=0,n_omega_exa
        w=w_omega(i_omega)
        do i_tau=0,n_tau/2
          c1=w*tau_mesh(i_tau)
          tau_from_omega(i_omega,i_tau,1)=cos(c1)
          tau_from_omega(i_omega,i_tau,2)=sin(c1)
        enddo
      enddo
c ---------- Geometric mesh spline contribution --------------------
      allocate(tmp(n_omega_exa-1:n_omega_exa+n_omega_geom+2,0:n_tau/2,
     &             2))
      call tau_from_freq_spl(n_tau/2,tau_mesh,n_omega_geom+3,
     &                       w_omega(n_omega_exa-1),tmp,0)
      do i_omega=n_omega_exa-1,n_omega_exa+n_omega_geom+2
        tau_from_omega(i_omega,:,1)=tau_from_omega(i_omega,:,1)
     &                             +tmp(i_omega,:,1)
        tau_from_omega(i_omega,:,2)=tau_from_omega(i_omega,:,2)
     &                             +tmp(i_omega,:,2)
      enddo
      deallocate(tmp)
c ---------- Semi-asymptotic spline contribution --------------------
      allocate(tmp(n_omega_exa+n_omega_geom:n_omega,0:n_tau/2,2))
      call tau_from_freq_spl(n_tau/2,tau_mesh,n_omega_asy+1,
     &                       w_omega(n_omega_exa+n_omega_geom),tmp,2)
      do i_omega=n_omega_exa+n_omega_geom,n_omega
        tau_from_omega(i_omega,:,1)=tau_from_omega(i_omega,:,1)
     &                             +tmp(i_omega,:,1)
        tau_from_omega(i_omega,:,2)=tau_from_omega(i_omega,:,2)
     &                             +tmp(i_omega,:,2)
      enddo
      deallocate(tmp)
c ---- Here are Omega-Sums for Internal Energy evaluation ---------
      allocate(et_sum_sc_gx(0:n_omega))
      allocate(et_sum_sc_gc(0:n_omega))
	et_sum_sc_gx=tau_from_omega(:,0,1)
      et_sum_sc_gc=et_sum_sc_gx  !! they differ only in asymptotic
c ---------- Asymptotic contribution --------------------------------
      tau_from_omega(n_omega,:,1)=tau_from_omega(n_omega,:,1)
     &                           +sum_omega_asy(:,2,1)*omega_max**4
      tau_from_omega(n_omega,:,2)=tau_from_omega(n_omega,:,2)
     &                           +sum_omega_asy(:,2,2)*omega_max**3
c ---------- Normalization -----------------------------------------
      tau_from_omega=tau_from_omega*4.d0/betta_t
c ------- Asymptotic contribution for Etotal = Sigma_C*G_X -------------
      et_sum_sc_gx(n_omega)=et_sum_sc_gx(n_omega)
     &                     +sum_omega_asy(0,1,1)*omega_max**2
c ------- Asymptotic contribution for Etotal = Sigma_C*G_C -------------
      et_sum_sc_gc(n_omega)=et_sum_sc_gc(n_omega)
     &                     +sum_omega_asy(0,2,1)*omega_max**4
c ---------- Normalization -----------------------------------------
      et_sum_sc_gc=et_sum_sc_gc/betta_t
      et_sum_sc_gx=et_sum_sc_gx/betta_t
      end
      
      
      
      
      subroutine transf_dtau_from_omega_spl
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: i_omega
      real*8 :: w
      real*8, allocatable :: tmp(:,:)
c ----------------------------------------------------------------
      dtau_from_omega=0.d0
c ---------- Exact contribution ----------------------------------  
      do i_omega=0,n_omega_exa
        w=w_omega(i_omega)
        dtau_from_omega(i_omega,1)=-(-1)**i_omega*w
        dtau_from_omega(i_omega,2)=w
      enddo
c ---------- Geometric mesh spline contribution --------------------
      allocate(tmp(n_omega_exa-1:n_omega_exa+n_omega_geom+2,2))
      call dtau_from_freq_spl(n_omega_geom+3,w_omega(n_omega_exa-1),tmp,
     &                        0)
      do i_omega=n_omega_exa-1,n_omega_exa+n_omega_geom+2
        dtau_from_omega(i_omega,1)=dtau_from_omega(i_omega,1)
     &                            +tmp(i_omega,1)
        dtau_from_omega(i_omega,2)=dtau_from_omega(i_omega,2)
     &                            +tmp(i_omega,2)
      enddo
      deallocate(tmp)
c ---------- Semi-asymptotic spline contribution --------------------
      allocate(tmp(n_omega_exa+n_omega_geom:n_omega,2))
      call dtau_from_freq_spl(n_omega_asy+1,
     &                        w_omega(n_omega_exa+n_omega_geom),tmp,2)
      do i_omega=n_omega_exa+n_omega_geom,n_omega
        dtau_from_omega(i_omega,1)=dtau_from_omega(i_omega,1)
     &                            +tmp(i_omega,1)
        dtau_from_omega(i_omega,2)=dtau_from_omega(i_omega,2)
     &                            +tmp(i_omega,2)
      enddo
      deallocate(tmp)
c ---------- Asymptotic contribution --------------------------------
      dtau_from_omega(n_omega,1)=dtau_from_omega(n_omega,1)
     &                          -sum_omega_asy(n_tau/2,2,2)*omega_max**4
      dtau_from_omega(n_omega,2)=dtau_from_omega(n_omega,2)
     &                          +sum_omega_asy(0,1,1)*omega_max**3
c ---------- Normalization -----------------------------------------
      dtau_from_omega=dtau_from_omega*4.d0/betta_t
      end
      
      
      subroutine transf_tau_from_omega_spl_3
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: i_omega,i_tau
      real*8 :: c1,w
      real*8, allocatable :: tmp(:,:,:)
c ----------------------------------------------------------------
      tau_from_omega=0.d0
c ---------- Exact contribution ----------------------------------  
      do i_omega=0,n_omega_exa
        w=w_omega(i_omega)
        do i_tau=0,n_tau/2
          c1=w*tau_mesh(i_tau)
          tau_from_omega(i_omega,i_tau,1)=cos(c1)
          tau_from_omega(i_omega,i_tau,2)=sin(c1)
        enddo
      enddo
c ---------- Geometric mesh spline contribution --------------------
      allocate(tmp(n_omega_exa-1:n_omega_exa+n_omega_geom+2,0:n_tau/2,
     &             2))
      call tau_from_freq_spl(n_tau/2,tau_mesh,n_omega_geom+3,
     &                       w_omega(n_omega_exa-1),tmp,0)
      do i_omega=n_omega_exa-1,n_omega_exa+n_omega_geom+2
        tau_from_omega(i_omega,:,1)=tau_from_omega(i_omega,:,1)
     &                             +tmp(i_omega,:,1)
        tau_from_omega(i_omega,:,2)=tau_from_omega(i_omega,:,2)
     &                             +tmp(i_omega,:,2)
      enddo
      deallocate(tmp)
c ---------- Semi-asymptotic spline contribution --------------------
      allocate(tmp(n_omega_exa+n_omega_geom:n_omega,0:n_tau/2,2))
      call tau_from_freq_spl(n_tau/2,tau_mesh,n_omega_asy+1,
     &                       w_omega(n_omega_exa+n_omega_geom),tmp,3)
      do i_omega=n_omega_exa+n_omega_geom,n_omega
        tau_from_omega(i_omega,:,1)=tau_from_omega(i_omega,:,1)
     &                             +tmp(i_omega,:,1)
        tau_from_omega(i_omega,:,2)=tau_from_omega(i_omega,:,2)
     &                             +tmp(i_omega,:,2)
      enddo
      deallocate(tmp)
c ---- Here are Omega-Sums for Internal Energy evaluation ---------
      allocate(et_sum_sc_gx(0:n_omega))
      allocate(et_sum_sc_gc(0:n_omega))
	et_sum_sc_gx=tau_from_omega(:,0,1)
      et_sum_sc_gc=et_sum_sc_gx  !! they differ only in asymptotic
c ---------- Asymptotic contribution --------------------------------
      tau_from_omega(n_omega,:,1)=tau_from_omega(n_omega,:,1)
     &                           +sum_omega_asy(:,2,1)*omega_max**4
      tau_from_omega(n_omega,:,2)=tau_from_omega(n_omega,:,2)
     &                           +sum_omega_asy(:,2,2)*omega_max**3
c ---------- Normalization -----------------------------------------
      tau_from_omega=tau_from_omega*4.d0/betta_t
c ------- Asymptotic contribution for Etotal = Sigma_C*G_X -------------
      et_sum_sc_gx(n_omega)=et_sum_sc_gx(n_omega)
     &                     +sum_omega_asy(0,1,1)*omega_max**2
c ------- Asymptotic contribution for Etotal = Sigma_C*G_C -------------
      et_sum_sc_gc(n_omega)=et_sum_sc_gc(n_omega)
     &                     +sum_omega_asy(0,2,1)*omega_max**4
c ---------- Normalization -----------------------------------------
      et_sum_sc_gc=et_sum_sc_gc/betta_t
      et_sum_sc_gx=et_sum_sc_gx/betta_t
      end
      
      
      
      
      subroutine transf_dtau_from_omega_spl_3
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: i_omega
      real*8 :: w
      real*8, allocatable :: tmp(:,:)
c ----------------------------------------------------------------
      dtau_from_omega=0.d0
c ---------- Exact contribution ----------------------------------  
      do i_omega=0,n_omega_exa
        w=w_omega(i_omega)
        dtau_from_omega(i_omega,1)=-(-1)**i_omega*w
        dtau_from_omega(i_omega,2)=w
      enddo
c ---------- Geometric mesh spline contribution --------------------
      allocate(tmp(n_omega_exa-1:n_omega_exa+n_omega_geom+2,2))
      call dtau_from_freq_spl(n_omega_geom+3,w_omega(n_omega_exa-1),tmp,
     &                        0)
      do i_omega=n_omega_exa-1,n_omega_exa+n_omega_geom+2
        dtau_from_omega(i_omega,1)=dtau_from_omega(i_omega,1)
     &                            +tmp(i_omega,1)
        dtau_from_omega(i_omega,2)=dtau_from_omega(i_omega,2)
     &                            +tmp(i_omega,2)
      enddo
      deallocate(tmp)
c ---------- Semi-asymptotic spline contribution --------------------
      allocate(tmp(n_omega_exa+n_omega_geom:n_omega,2))
      call dtau_from_freq_spl(n_omega_asy+1,
     &                        w_omega(n_omega_exa+n_omega_geom),tmp,3)
      do i_omega=n_omega_exa+n_omega_geom,n_omega
        dtau_from_omega(i_omega,1)=dtau_from_omega(i_omega,1)
     &                            +tmp(i_omega,1)
        dtau_from_omega(i_omega,2)=dtau_from_omega(i_omega,2)
     &                            +tmp(i_omega,2)
      enddo
      deallocate(tmp)
c ---------- Asymptotic contribution --------------------------------
      dtau_from_omega(n_omega,1)=dtau_from_omega(n_omega,1)
     &                          -sum_omega_asy(n_tau/2,2,2)*omega_max**4
      dtau_from_omega(n_omega,2)=dtau_from_omega(n_omega,2)
     &                          +sum_omega_asy(0,1,1)*omega_max**3
c ---------- Normalization -----------------------------------------
      dtau_from_omega=dtau_from_omega*4.d0/betta_t
      end