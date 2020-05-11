      subroutine vertex_tau_to_omega_uniform(lambda_tau,lambda_dyn,isp,
     &                                       jsp)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: isp,jsp
      complex*16, intent(in) :: lambda_tau(nbndf_bnd,nbndf_bnd,nqdiv_c,
     &                                     2,0:n_tau/2)
      complex*16, intent(out) :: lambda_dyn(nbndf_bnd,nbndf_bnd,
     &                                      0:n_omega,nqdiv_c)
      integer :: k0,nk,k,kf0,it,i_tau
      complex*16, allocatable :: vxt(:,:,:,:),lambda(:,:,:,:)
c -------- Transform Vrt(tau;nu) ---- > Vrt(omega;nu) ------------
	do k=1,nqdiv_c
	  k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	  nk=n_low_bnd(kf0,isp)
	  allocate(vxt(nk,nk,2,0:n_tau/2))
	  allocate(lambda(nk,nk,0:n_omega,2))
        do i_tau=0,n_tau/2
          do it=1,2
            vxt(:,:,it,i_tau)=lambda_tau(1:nk,1:nk,k,it,i_tau)
          enddo
        enddo
        call tau_to_omega_kk(lambda,vxt,nk*nk)
        lambda_dyn(1:nk,1:nk,:,k)=lambda(:,:,:,2) 
        deallocate(vxt,lambda)
      enddo   !! over k
      end
