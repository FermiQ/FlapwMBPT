      subroutine g_p_bnd_mm(g_tau,g_omega)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  real*8, intent(out) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                             npnt_c,nspin),
     &                       g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,
     &                               npnt_c,nspin)
      integer :: k,ispin,ib,ind_k,ind_tau,it,jbnd,jb,ibnd,kc,n
      real*8, allocatable :: grt(:,:,:,:)
      g_tau=0.d0
      g_omega=0.d0
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          kc=k_c_from_a(k)
          if(kc==0) cycle
          n=n_low_bnd(k,ispin)
          allocate(grt(n,n,2,ndim3_tau))
          do ind_tau=1,ndim3_tau
            do it=1,2
	          do jbnd=1,n
	            jb=ind_bands_bnd(jbnd,k,ispin)
	            do ibnd=1,n
	              ib=ind_bands_bnd(ibnd,k,ispin)
                  grt(ibnd,jbnd,it,ind_tau)=
     &	              g_full(ib,jb,it,ind_tau,ind_k,ispin)
                enddo
              enddo
            enddo
          enddo   !! over ind_tau
          call tau_to_omega_baa(g_omega(1,1,1,1,kc,ispin),nbndf_bnd,grt,
     &                          n,n)
          do ind_tau=1,ndim3_tau
            g_tau(1:n,1:n,:,ind_tau,kc,ispin)=grt(:,:,:,ind_tau)
          enddo
          deallocate(grt)
	    enddo  ! over ind_k
	  enddo
	  if(nproc_k/=1) then
	    call dgop(g_tau,2*nbndf_bnd**2*npnt_c*ndim3_tau*nspin,'  +',
     &            comm_k)
        call dgop(g_omega,2*nbndf_bnd**2*npnt_c*ndim3_omega*nspin,'  +',
     &            comm_k)
      endif
      end
