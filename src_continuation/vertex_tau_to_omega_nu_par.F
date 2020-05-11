      subroutine vertex_tau_to_omega_nu_par(lambda_tau,lambda_dyn,i_nu,
     &                                      iq,nom,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_nu,iq,nom,ispin
      complex*16, intent(in) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                     ndimc_kk(me_k+1),2,2,
     &                                     ndim3_tau)
      complex*16, intent(out) :: lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,
     &                                      ndimc_kk(me_k+1))
      integer :: k0,nk,kq,kq0,k,nkq,il,kf0,kfq0,it,i_tau,ind_tau,ind_k
      real*8 :: v(3)
      complex*16, allocatable :: vxt(:,:,:,:,:),lambda(:,:,:,:)
c -------- Transform Vrt(tau;nu) ---- > Vrt(omega;nu) ------------
      do ind_k=1,ndimc_kk(me_k+1)
        k=nc_mpi_kk(me_k+1)+ind_k
        k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
        nk=n_low_bnd(kf0,ispin)
        v=pnt_c(:,k)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)     !! for K-Q
        kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
        nkq=n_low_bnd(kfq0,ispin)
        allocate(vxt(nk,nkq,2,ndim3_tau,2))
        allocate(lambda(nk,nkq,nom,2))
        do ind_tau=1,ndim3_tau
          do il=1,2
            do it=1,2
              vxt(:,:,it,ind_tau,il)=
     &          lambda_tau(1:nk,1:nkq,ind_k,il,it,ind_tau)
            enddo
          enddo
        enddo
        do il=1,2
          call from_tau_nu_to_omega_nu(nk*nkq,i_nu,nom,lambda(1,1,1,il),
     &                                 vxt,il)
          lambda_dyn(1:nk,1:nkq,1:nom,il,ind_k)=lambda(:,:,1:nom,il)
        enddo
        deallocate(vxt,lambda)
      enddo   !! over ind_k
      end
