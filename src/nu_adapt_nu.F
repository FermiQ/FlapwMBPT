      subroutine nu_adapt_nu
      use heg_mod
      use manager_mod
      use parallel_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: k,i,i_nu,i_tau,ind_tau,it,ind1_nu,i1_nu,ind_nu,nas
      real*8 :: hh,pib
      real*8, allocatable :: mesh(:)
      complex*16, allocatable :: tm(:,:)
      pib=pi/betta_t
      hh=pib+pib
      nas=n_nu_geom+n_nu_asy
      allocate(mesh(0:nas))
      allocate(num_nu2_adapt(0:n_nu))
      num_nu2_adapt=0
      n_nu2_max=2*(n_nu_geom+n_nu_asy+n_nu_exa)+2
      allocate(w_nu_adapt_nu(n_nu2_max,0:n_nu))
      allocate(tau_nu_from_nu_nu(n_nu2_max,2,ndim3_tau,
     &	                       ndim_nu(me_k+1)))
      tau_nu_from_nu_nu=(0.d0,0.d0)
      w_nu_adapt_nu=0.d0
      do ind1_nu=1,ndim_nu(me_k+1)
        i1_nu=n_mpi_nu(me_k+1)+ind1_nu-1
        call nu2_mesh_spl(w_nu(i1_nu),w_nu_adapt_nu(1,i1_nu),
     &                    tau_nu_from_nu_nu(1,1,1,ind1_nu),
     &                    num_nu2_adapt(i1_nu))
      enddo   !! over ind1_nu
      if(nproc_k/=1) then
        call igop(num_nu2_adapt,n_nu+1,'  +',comm_k)
        call dgop(w_nu_adapt_nu,n_nu2_max*(n_nu+1),'  +',comm_k)
      endif
      tau_nu_from_nu_nu=tau_nu_from_nu_nu/betta_t
      deallocate(mesh)
c ---------------------------------------------------------------------
      allocate(nu_nu_from_tau_nu(0:n_tau/2,n_nu2_max,2,2,ndim3_nu))
      nu_nu_from_tau_nu=0.d0
      allocate(mesh(n_nu2_max))
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        do i=1,num_nu2_adapt(i_nu)
          mesh(i)=w_nu_adapt_nu(i,i_nu)
        enddo
        call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),
     &                               num_nu2_adapt(i_nu)-1,mesh,
     &                        nu_nu_from_tau_nu(0,1,1,1,ind_nu),
     &                               n_nu2_max-1,0)
        do i=1,num_nu2_adapt(i_nu)
          mesh(i)=w_nu_adapt_nu(i,i_nu)-w_nu(i_nu)
        enddo
        call transf_bos_from_tau_spl(n_tau/2,tau_mesh(0),
     &                               num_nu2_adapt(i_nu)-1,mesh,
     &                        nu_nu_from_tau_nu(0,1,1,2,ind_nu),
     &                               n_nu2_max-1,0)
      enddo  !! over ind_nu
      deallocate(mesh)
c -----------------------------------------------------------	      
      allocate(tau_nu_from_nu_nu_1(n_nu2_max,0:n_tau,ndim3_nu))
      allocate(tm(0:n_tau,0:n_nu))
      do k=1,n_nu2_max
        tm=(0.d0,0.d0)
        do ind_nu=1,ndim_nu(me_k+1)
          i_nu=n_mpi_nu(me_k+1)+ind_nu-1
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=ndim3_tau*me_t+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
              tm(i_tau,i_nu)=tau_nu_from_nu_nu(k,it,ind_tau,ind_nu)
            enddo
          enddo
        enddo
        if(nproc_t/=1) call dgop(tm,2*(n_tau+1)*(n_nu+1),'  +',
     &	                         comm_t)
        if(nproc_k/=1) then
          call dgop(tm,2*(n_tau+1)*(n_nu+1),'  +',comm_k)
        endif
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          do i_tau=0,n_tau
            tau_nu_from_nu_nu_1(k,i_tau,ind_nu)=conjg(tm(i_tau,i_nu))
          enddo
        enddo
      enddo
      deallocate(tm)
      end
