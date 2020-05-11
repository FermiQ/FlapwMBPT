      subroutine omega_adapt_omega
      use heg_mod
      use manager_mod
      use parallel_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i_omega,k,i_tau,ind_omega,ind_tau,it,nas,i
      real*8 :: hh,pib
      real*8, allocatable :: mesh(:)
      complex*16, allocatable :: tm(:,:)
      pib=pi/betta_t
      hh=pib+pib
      nas=n_omega_geom+n_omega_asy
      allocate(num_omega1_omega(0:n_omega))
      num_omega1_omega=0
      n_omega1_max=2*(n_omega_asy+n_omega_geom+n_omega_exa)+2
      allocate(w_om_adapt_omega(n_omega1_max,0:n_omega))
      allocate(tau_from_omega_omega(n_omega1_max,2,ndim3_tau,0:n_omega))
      tau_from_omega_omega=(0.d0,0.d0)
      do i_omega=0,n_omega
        call omega_adapt_nu_spl(w_omega(i_omega)-pib,
     &                          w_om_adapt_omega(1,i_omega),
     &                          tau_from_omega_omega(1,1,1,i_omega),
     &                          num_omega1_omega(i_omega))
      enddo   !! over i_nu
      tau_from_omega_omega=tau_from_omega_omega/betta_t
c --------------------------------------------------------------
      allocate(tau_om_from_omega_om(n_omega1_max,0:n_tau,ndim3_omega))
      allocate(tm(0:n_tau,0:n_omega))
      do k=1,n_omega1_max
        tm=(0.d0,0.d0)
        do i_omega=0,n_omega
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=me_t*ndim3_tau+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
              tm(i_tau,i_omega)=
     &          tau_from_omega_omega(k,it,ind_tau,i_omega)
            enddo
          enddo
        enddo
        if(nproc_t/=1) call dgop(tm,2*(n_tau+1)*(n_omega+1),'  +',
     &	                         comm_t)
        do ind_omega=1,ndim3_omega
          i_omega=me_t*ndim3_omega+ind_omega-1
          do i_tau=0,n_tau
            tau_om_from_omega_om(k,i_tau,ind_omega)=tm(i_tau,i_omega)
          enddo
        enddo
      enddo
      deallocate(tm)
c ---------------------------------------------------------------------
      allocate(omega_om_from_tau_om(0:n_tau/2,n_omega1_max,4,2,
     &                              ndim3_omega))
      omega_om_from_tau_om=0.d0
      allocate(mesh(n_omega1_max))
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        do i=1,num_omega1_omega(i_omega)
          mesh(i)=w_om_adapt_omega(i,i_omega)+pib
        enddo
        call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),
     &                               num_omega1_omega(i_omega)-1,mesh,
     &                        omega_om_from_tau_om(0,1,1,1,ind_omega),
     &                               n_omega1_max-1,0)
        do i=1,num_omega1_omega(i_omega)
          mesh(i)=w_om_adapt_omega(i,i_omega)-w_omega(i_omega)+pib+pib
        enddo
        call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),
     &                               num_omega1_omega(i_omega)-1,mesh,
     &                        omega_om_from_tau_om(0,1,3,1,ind_omega),
     &                               n_omega1_max-1,0)
        do i=1,num_omega1_omega(i_omega)
          mesh(i)=w_om_adapt_omega(i,i_omega)-pib
        enddo
        call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),
     &                               num_omega1_omega(i_omega)-1,mesh,
     &                        omega_om_from_tau_om(0,1,1,2,ind_omega),
     &                               n_omega1_max-1,0)
        do i=1,num_omega1_omega(i_omega)
          mesh(i)=-w_om_adapt_omega(i,i_omega)+w_omega(i_omega)
        enddo
        call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),
     &                               num_omega1_omega(i_omega)-1,mesh,
     &                        omega_om_from_tau_om(0,1,3,2,ind_omega),
     &                               n_omega1_max-1,0)
      enddo  !! over ind_omega
      deallocate(mesh)
      end
