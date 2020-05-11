      subroutine vertex_tau_to_omega_loc(lambda_tau,lambda_dyn,ind_nu,
     &                                   nom,isp,jsp,ll,n)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,nom,isp,jsp,ll,n
      complex*16, intent(in) :: lambda_tau(n,n,ll,2,0:n_tau/2)
      complex*16, intent(out) :: lambda_dyn(n,n,nom,2)
      integer :: il,it,i_tau
      complex*16, allocatable :: vxt(:,:,:,:,:),lambda(:,:,:,:)
c -------- Transform Vrt(tau;nu) ---- > Vrt(omega;nu) ------------
      allocate(vxt(n,n,2,0:n_tau/2,ll))
      allocate(lambda(n,n,nom,2))
      do i_tau=0,n_tau/2
        do il=1,ll
          do it=1,2
            vxt(:,:,it,i_tau,il)=lambda_tau(:,:,il,it,i_tau)
          enddo
        enddo
      enddo
      if(ll==1) then
        call tau_to_omega_kk(lambda,vxt,n*n)
        lambda_dyn(:,:,1:nom,1)=lambda(:,:,1:nom,1)
        lambda_dyn(:,:,1:nom,2)=lambda(:,:,1:nom,2)
      else if(ll==2) then
        do il=1,2
          call from_tau_nu_to_omega_nu_kk(n*n,ind_nu,nom,
     &                                    lambda(1,1,1,il),vxt,il)
          lambda_dyn(:,:,1:nom,il)=lambda(:,:,1:nom,il)
        enddo
      endif
      deallocate(vxt,lambda)
      end
