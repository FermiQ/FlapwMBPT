      subroutine transf_tau_from_nu_spl
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: i_nu,i_tau
      real*8 :: c1,w
      real*8, allocatable :: tmp(:,:,:)
c ----------------------------------------------------------------
      tau_from_nu=0.d0
c ---------- Exact contribution ----------------------------------      
      tau_from_nu(0,:,1)=0.5d0
      do i_nu=1,n_nu_exa
        w=w_nu(i_nu)
        do i_tau=0,n_tau/2
          c1=w*tau_mesh(i_tau)
          tau_from_nu(i_nu,i_tau,1)=cos(c1)
          tau_from_nu(i_nu,i_tau,2)=sin(c1)
        enddo
      enddo
c ---------- Geometric mesh spline contribution --------------------
      allocate(tmp(n_nu_exa-1:n_nu_exa+n_nu_geom+2,0:n_tau/2,2))
      call tau_from_freq_spl(n_tau/2,tau_mesh,n_nu_geom+3,
     &                       w_nu(n_nu_exa-1),tmp,0)
      do i_nu=n_nu_exa-1,n_nu_exa+n_nu_geom+2
        tau_from_nu(i_nu,:,1)=tau_from_nu(i_nu,:,1)+tmp(i_nu,:,1)
        tau_from_nu(i_nu,:,2)=tau_from_nu(i_nu,:,2)+tmp(i_nu,:,2)
      enddo
      deallocate(tmp)
c ---------- Semi-asymptotic spline contribution --------------------
      allocate(tmp(n_nu_exa+n_nu_geom:n_nu,0:n_tau/2,2))
      call tau_from_freq_spl(n_tau/2,tau_mesh,n_nu_asy+1,
     &                       w_nu(n_nu_exa+n_nu_geom),tmp,2)
      do i_nu=n_nu_exa+n_nu_geom,n_nu
        tau_from_nu(i_nu,:,1)=tau_from_nu(i_nu,:,1)+tmp(i_nu,:,1)
        tau_from_nu(i_nu,:,2)=tau_from_nu(i_nu,:,2)+tmp(i_nu,:,2)
      enddo
      deallocate(tmp)
c ---------- Asymptotic contribution --------------------------------
      tau_from_nu(n_nu,:,1)=tau_from_nu(n_nu,:,1)
     &                     +sum_nu_asy(:,1,1)*nu_max**2
      tau_from_nu(n_nu,:,2)=tau_from_nu(n_nu,:,2)
     &                     +sum_nu_asy(:,2,2)*nu_max**3
c ---- Here are Nu-Sums for Internal Energy evaluation ---------
      allocate(et_sum_nu(0:n_nu))
      et_sum_nu=tau_from_nu(:,0,1)
c ------- Asymptotic contribution for Etotal = P*W --------------------
      et_sum_nu(n_nu)=et_sum_nu(n_nu)+sum_nu_asy(0,2,1)*nu_max**4
c ---------- Normalization -----------------------------------------
      et_sum_nu=et_sum_nu/betta_t
      tau_from_nu=tau_from_nu*4.d0/betta_t
      end