      subroutine vertex_gt_loc(it,i0_tau,ll,q_pw,fif,lambda_tau,n,np,
     &                         gtau)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,i0_tau,ll,n,np
      real*8, intent(in) :: fif(nrel,n,n,np,nspin),
     &                      gtau(nrel,n,n,0:n_tau,nspin)
      complex*16, intent(in) :: q_pw(np,np,ll,2)
      complex*16, intent(inout) :: lambda_tau(n,n,ll,2,0:n_tau/2,
     &                                        nspin)
      integer :: i_tau,ispin,il
      complex*16, allocatable :: vx(:,:,:,:)
      i_tau=i0_tau
      if(it==2) i_tau=n_tau-i0_tau
      allocate(vx(n,n,ll,nspin))
      vx=(0.d0,0.d0)
      if(irel/=2) then
        do ispin=1,nspin
          call vertex_loc_gt_r(i_tau,ll,q_pw(1,1,1,it),
     &                          gtau(1,1,1,0,ispin),
     &                          vx(1,1,1,ispin),fif(1,1,1,1,ispin),n,np)
        enddo
      else
        call vertex_loc_gt_c(i_tau,ll,q_pw(1,1,1,it),gtau,vx,fif,n,np)
      endif
      do ispin=1,nspin
        do il=1,ll
          lambda_tau(:,:,il,it,i0_tau,ispin)=
     &      lambda_tau(:,:,il,it,i0_tau,ispin)+vx(:,:,il,ispin)
        enddo
      enddo
      deallocate(vx)
      end
