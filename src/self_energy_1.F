      subroutine self_energy_1
	  use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
	  use units_mod
	  use vertex_mod
	  implicit none
	  integer :: ispin,ind_k,k,i,it,ind_omega,n,j,j0,i0,kc,i_omega
	  real*8, allocatable :: sigc(:,:,:,:),sd(:,:,:,:,:),zd(:,:,:)
      complex*16, allocatable :: sigc0(:,:),zren(:,:),tmp(:,:)
	  if(ubi==' qp'.or.ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
	    if(iter_psi/=0.or.iter_bsp/=0) then
	      if(allocated(z_wi)) deallocate(z_wi)
	      allocate(z_wi(nbndf_bnd,nbndf_bnd,npnt_c,nspin))
	      z_wi=(0.d0,0.d0)
          if(allocated(sig_c_wi)) deallocate(sig_c_wi)
          allocate(sig_c_wi(nbndf_bnd,nbndf_bnd,npnt_c,0:1,nspin))
          sig_c_wi=(0.d0,0.d0)
	    endif
	    allocate(sd(nbndf,2,ndim3_omega,ndim3_k(me_k+1),nspin))
        allocate(zd(nbndf,ndim3_k(me_k+1),nspin))
c ------- We prepare Sigma and its derivative at Omega=0 -----------
	    do ispin=1,nspin
	      do ind_k=1,ndim3_k(me_k+1)
	        k=n3_mpi_k(me_k+1)+ind_k
	        n=n_bnd(k,ispin)
            allocate(sigc(n,n,2,ndim3_omega))
            call tau_to_omega_baa(sigc,n,sig_c_tau(1,1,1,1,ind_k,ispin),
     &                            nbndf,n)
            do ind_omega=1,ndim3_omega
              do it=1,2
                do i=1,n
                  sd(i,it,ind_omega,ind_k,ispin)=sigc(i,i,it,ind_omega)
                enddo
              enddo
            enddo
            allocate(sigc0(n,n))
            allocate(zren(n,n))
	        call sig_omega0_solid(n,sigc,ind_k,ispin,sigc0,zren)
	        deallocate(sigc0)
	        if(iter_psi/=0.or.iter_bsp/=0) then
	          kc=k_c_from_a(k)
	          if(kc/=0) then
                do j=1,n_low_bnd(k,ispin)
                  j0=ind_bands_bnd(j,k,ispin)
                  do i=1,n_low_bnd(k,ispin)
                    i0=ind_bands_bnd(i,k,ispin)
                    z_wi(i,j,kc,ispin)=zren(i0,j0)
                    if(i==j) z_wi(i,j,kc,ispin)=z_wi(i,j,kc,ispin)-1.d0
                  enddo
                enddo
                allocate(tmp(n,n))
                do ind_omega=1,ndim3_omega
                  i_omega=me_t*ndim3_omega+ind_omega-1
                  if(i_omega<=1) then
                    call ferm_unpack_omega(tmp,sigc(1,1,1,ind_omega),n,
     &                                     n,n)
                    do j=1,n_low_bnd(k,ispin)
                      j0=ind_bands_bnd(j,k,ispin)
                      do i=1,n_low_bnd(k,ispin)
                        i0=ind_bands_bnd(i,k,ispin)
                        sig_c_wi(i,j,kc,i_omega,ispin)=tmp(i0,j0)
                      enddo
                    enddo
                  endif
                enddo
                deallocate(tmp)
              endif
            endif
            deallocate(sigc)
	        call invers_h(n,zren,n)
	        do i=1,n
	          zd(i,ind_k,ispin)=zren(i,i)
	        enddo
	        deallocate(zren)
	      enddo
	    enddo
c ------------------------------------------------------------------
        call output_z_factor(zd)
        call output_sigma_band(sd,'sigc_w','sigc',evolt/2)
        deallocate(sd,zd)
        if(iter_psi/=0.or.iter_bsp/=0) then
          if(nproc_k/=1) then
            call dgop(z_wi,2*nbndf_bnd**2*npnt_c*nspin,'  +',comm_k)
            call dgop(sig_c_wi,4*nbndf_bnd**2*npnt_c*nspin,'  +',comm_k)
          endif
          if(nproc_t/=1) then
            call dgop(sig_c_wi,4*nbndf_bnd**2*npnt_c*nspin,'  +',comm_t)
          endif
        endif
      endif
c ------------------------------------------------------------------
      call timel('**** SELF_ENERGY_1 finished ********')
      end
