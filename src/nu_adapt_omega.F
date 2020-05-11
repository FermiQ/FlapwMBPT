      subroutine nu_adapt_omega
      use heg_mod
      use manager_mod
      use parallel_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: k,i,i_omega,i_tau,ind_tau,it,ind_omega,nas
      real*8 :: hh,pib
      real*8, allocatable :: mesh(:)
      complex*16, allocatable :: tm(:,:)
      pib=pi/betta_t
      hh=pib+pib
      nas=n_nu_geom+n_nu_asy
      allocate(mesh(0:nas))
      allocate(num_nu_omega(0:n_omega))
      num_nu_omega=0
      n_nu2_max=2*(n_nu_geom+n_nu_asy+n_nu_exa)+2
      allocate(w_nu_adapt_omega(n_nu2_max,0:n_omega))
      allocate(tau_om_from_nu_om(n_nu2_max,2,ndim3_tau,
     &	                       ndim_omega(me_k+1)))
      tau_om_from_nu_om=(0.d0,0.d0)
      w_nu_adapt_omega=0.d0
      do ind_omega=1,ndim_omega(me_k+1)
        i_omega=n_mpi_omega(me_k+1)+ind_omega-1
        call nu2_mesh_spl(w_omega(i_omega)-pib,
     &                    w_nu_adapt_omega(1,i_omega),
     &                    tau_om_from_nu_om(1,1,1,ind_omega),
     &                    num_nu_omega(i_omega))
      enddo   !! over ind_omega
      if(nproc_k/=1) then
        call igop(num_nu_omega,n_omega+1,'  +',comm_k)
        call dgop(w_nu_adapt_omega,n_nu2_max*(n_omega+1),'  +',comm_k)
      endif
      tau_om_from_nu_om=tau_om_from_nu_om/betta_t
      deallocate(mesh)
c ---------------------------------------------------------------------
      allocate(nu_om_from_tau_om(0:n_tau/2,n_nu2_max,2,2,ndim3_omega))
      nu_om_from_tau_om=0.d0
      allocate(mesh(n_nu2_max))
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        do i=1,num_nu_omega(i_omega)
          mesh(i)=w_nu_adapt_omega(i,i_omega)
        enddo
        call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),
     &                               num_nu_omega(i_omega)-1,mesh,
     &                            nu_om_from_tau_om(0,1,1,1,ind_omega),
     &                               n_nu2_max-1,0)
        do i=1,num_nu_omega(i_omega)
          mesh(i)=w_nu_adapt_omega(i,i_omega)-w_omega(i_omega)+pib
        enddo
        call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),
     &                               num_nu_omega(i_omega)-1,mesh,
     &                        nu_om_from_tau_om(0,1,1,2,ind_omega),
     &                               n_nu2_max-1,0)
      enddo  !! over ind_omega
      deallocate(mesh)
c -----------------------------------------------------------	      
      allocate(tau_om_from_nu_om_1(n_nu2_max,0:n_tau,ndim3_omega))
      allocate(tm(0:n_tau,0:n_omega))
      do k=1,n_nu2_max
        tm=(0.d0,0.d0)
        do ind_omega=1,ndim_omega(me_k+1)
          i_omega=n_mpi_omega(me_k+1)+ind_omega-1
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=ndim3_tau*me_t+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
              tm(i_tau,i_omega)=
     &          tau_om_from_nu_om(k,it,ind_tau,ind_omega)
            enddo
          enddo
        enddo
        if(nproc_t/=1) call dgop(tm,2*(n_tau+1)*(n_omega+1),'  +',
     &	                         comm_t)
        if(nproc_k/=1) then
          call dgop(tm,2*(n_tau+1)*(n_omega+1),'  +',comm_k)
        endif
        do ind_omega=1,ndim3_omega
          i_omega=me_t*ndim3_omega+ind_omega-1
          do i_tau=0,n_tau
            tau_om_from_nu_om_1(k,i_tau,ind_omega)=
     &        conjg(tm(i_tau,i_omega))
          enddo
        enddo
      enddo
      deallocate(tm)
      end
