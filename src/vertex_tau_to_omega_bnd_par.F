      subroutine vertex_tau_to_omega_bnd_par(lambda_tau,lambda_dyn,
     &                                       ind_nu,iq,nom,isp,jsp,ll)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,iq,nom,isp,jsp,ll
      complex*16, intent(in) :: lambda_tau(nbndf_bnd,nbndf_bnd,nqdiv_c,
     &                                     ll,2,0:n_tau/2)
      complex*16, intent(out) :: lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,
     &                                      nqdiv_c)
      integer :: k0,nk,kq,kq0,k,nkq,il,kf0,kfq0,it,i_tau
      real*8 :: v(3)
      complex*16, allocatable :: vxt(:,:,:,:,:),lambda(:,:,:,:)
c -------- Transform Vrt(tau;nu) ---- > Vrt(omega;nu) ------------
	  do k=1,nqdiv_c
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,isp)
	    v=pnt_c(:,k)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K-Q
	    kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kfq0,jsp)
	    allocate(vxt(nk,nkq,2,0:n_tau/2,ll))
	    allocate(lambda(nk,nkq,nom,2))
        do i_tau=0,n_tau/2
          do il=1,ll
            do it=1,2
              vxt(:,:,it,i_tau,il)=
     &          lambda_tau(1:nk,1:nkq,k,il,it,i_tau)
            enddo
          enddo
        enddo
        if(ll==1) then
          call tau_to_omega_kk(lambda,vxt,nk*nkq)
          lambda_dyn(1:nk,1:nkq,1:nom,1,k)=lambda(:,:,1:nom,1)
          lambda_dyn(1:nk,1:nkq,1:nom,2,k)=lambda(:,:,1:nom,2) 
        else if(ll==2) then
          do il=1,2
            call from_tau_nu_to_omega_nu_kk(nk*nkq,ind_nu,nom,
     &                                     lambda(1,1,1,il),vxt,il)
            lambda_dyn(1:nk,1:nkq,1:nom,il,k)=lambda(:,:,1:nom,il)
          enddo
        endif
        deallocate(vxt,lambda)
      enddo   !! over k
      end
