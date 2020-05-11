      subroutine omega_adapt_nu
      use heg_mod
      use manager_mod
      use parallel_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i_omega,k,i_nu,i_tau,ind_nu,nom,ind_tau,it,nas
      real*8 :: hh,pib
      real*8, allocatable :: mesh(:),t(:,:)
      complex*16, allocatable :: tm(:,:)
      pib=pi/betta_t
      hh=pib+pib
      nas=n_omega_geom+n_omega_asy
      allocate(num_omega1_adapt(0:n_nu))
      num_omega1_adapt=0
      n_omega1_max=2*(n_omega_asy+n_omega_geom+n_omega_exa)+2
      allocate(w_om_adapt_nu(n_omega1_max,0:n_nu))
      allocate(tau_from_omega_nu(n_omega1_max,2,ndim3_tau,0:n_nu))
      tau_from_omega_nu=(0.d0,0.d0)
      do i_nu=0,n_nu
        call omega_adapt_nu_spl(w_nu(i_nu),w_om_adapt_nu(1,i_nu),
     &                          tau_from_omega_nu(1,1,1,i_nu),
     &                          num_omega1_adapt(i_nu))
      enddo   !! over i_nu
      tau_from_omega_nu=tau_from_omega_nu/betta_t
c --------------------------------------------------------------
      allocate(tau_nu_from_omega_nu(n_omega1_max,0:n_tau,ndim3_nu))
      allocate(tm(0:n_tau,0:n_nu))
      do k=1,n_omega1_max
        tm=(0.d0,0.d0)
        do i_nu=0,n_nu
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=me_t*ndim3_tau+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
              tm(i_tau,i_nu)=tau_from_omega_nu(k,it,ind_tau,i_nu)
            enddo
          enddo
        enddo
        if(nproc_t/=1) call dgop(tm,2*(n_tau+1)*(n_nu+1),'  +',
     &	                         comm_t)
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          do i_tau=0,n_tau
            tau_nu_from_omega_nu(k,i_tau,ind_nu)=tm(i_tau,i_nu)
          enddo
        enddo
      enddo
      deallocate(tm)
c ---   Smooth transform F(t;v) ---- > F(w;v) ----------------	
      allocate(omega_nu_from_tau_nu_ab(0:n_tau/2,n_omega1_max,4,
     &	                               ndim3_nu))
      omega_nu_from_tau_nu_ab=0.d0
      allocate(mesh(n_omega1_max))
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        nom=num_omega1_adapt(i_nu)
        do i_omega=1,nom
          mesh(i_omega)=w_om_adapt_nu(i_omega,i_nu)
        enddo
        call transf_ferm_from_tau_spl(n_tau/2,tau_mesh(0),nom-1,
     &                                mesh,
     &                        omega_nu_from_tau_nu_ab(0,1,1,ind_nu),
     &                                n_omega1_max-1,0)
        do i_omega=1,nom
          mesh(i_omega)=w_om_adapt_nu(i_omega,i_nu)-w_nu(i_nu)
        enddo
        call transf_ferm_from_tau_spl(n_tau/2,tau_mesh(0),nom-1,
     &                                mesh,
     &                        omega_nu_from_tau_nu_ab(0,1,3,ind_nu),
     &                                n_omega1_max-1,0)
      enddo  !! over ind_nu
      deallocate(mesh)
      allocate(omega_nu_from_tau_nu(ndim3_tau,n_omega1_max,4,0:n_nu))
      omega_nu_from_tau_nu=0.d0
      allocate(t(0:n_tau/2,0:n_nu))
      do it=1,4
        do k=1,n_omega1_max
          t=0.d0
          do ind_nu=1,ndim3_nu
            i_nu=me_t*ndim3_nu+ind_nu-1
            t(:,i_nu)=omega_nu_from_tau_nu_ab(:,k,it,ind_nu)
          enddo
          if(nproc_t/=1) call dgop(t,(n_tau/2+1)*(n_nu+1),'  +',comm_t)
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            do i_nu=0,n_nu
              omega_nu_from_tau_nu(ind_tau,k,it,i_nu)=t(i_tau,i_nu)
            enddo
          enddo
        enddo
      enddo
      deallocate(t)
      end
