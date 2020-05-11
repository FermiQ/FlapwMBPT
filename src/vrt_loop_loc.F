      subroutine vrt_loop_loc(spin_flips,i0_tau,isp,jsp,ll,kloc,
     &                        w_tau,fif,lambda_tau,n,np,gtau)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical, intent(in) :: spin_flips
      integer, intent(in) :: i0_tau,ll,isp,jsp,n,np
      real*8, intent(in) :: w_tau(np,np),fif(nrel,n,n,np,nspin),
     &                      gtau(nrel,n,n,0:n_tau)
      complex*16, intent(in) :: kloc(n,n,ll,2)
      complex*16, intent(inout) :: lambda_tau(n,n,ll,2)
      integer :: il,it,i_tau
      complex*16, allocatable :: vx(:,:)
      allocate(vx(n,n))
      do il=1,ll
        do it=1,2
          vx=(0.d0,0.d0)
          if(irel<2) call vertex_vk_0_r(w_tau,kloc(1,1,il,it),vx,
     &                                  fif(1,1,1,1,isp),n,np,np)
          if(irel==2) call vertex_vk_0_c(w_tau,kloc(1,1,il,it),vx,fif,n,
     &                                   np,np)
          lambda_tau(:,:,il,it)=vx
          if(spin_flips) then
            i_tau=i0_tau
            if(it==2) i_tau=n_tau-i_tau
            if(irel<2) call q_tau_loc_r(i_tau,ll,it,il,
     &                                  qt_loc(1,1,il,it,i0_tau),
     &                                  gtau,kloc,fif(1,1,1,1,isp),n,np)
            if(irel==2) call q_tau_loc_c(i_tau,ll,it,il,
     &                                   qt_loc(1,1,il,it,i0_tau),
     &                                   gtau,kloc,fif(1,1,1,1,isp),n,
     &                                   np)
          endif
        enddo
      enddo
      deallocate(vx)
      end
