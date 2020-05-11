      subroutine vertex_loc_loop(spin_flips,ind_nu,it_vrt,ll,s2,nom,
     &                           nomnu,k_pw,fif,vloc,wloc,w_tau,
     &                           lambda_stat,lambda_dyn,nsp,n,np,gtau)
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
      integer, intent(in) :: n,np,ind_nu,nom,nomnu,s2,it_vrt,ll,nsp
      real*8, intent(in) :: vloc(n_pbmtm_red,n_pbmtm_red),
     &                      wloc(n_pbmt_red,n_pbmt_red,0:n_nu),
     &                      w_tau(np,np,0:n_tau/2),
     &                      fif(nrel,n,n,np,nspin),
     &                      gtau(nrel,n,n,0:n_tau,nspin)
      complex*16, intent(in) :: k_pw(n,n,0:n_tau,ll,nsp)
      complex*16, intent(out) :: lambda_stat(n,n,nsp),
     &                           lambda_dyn(n,n,nom,2,nsp)
      integer :: i0_tau,ispin,it,i_nu,isp,jsp,il,i_tau
      complex*16, allocatable :: lambda_tau(:,:,:,:,:,:),
     &                           ktau0(:,:,:),kloc(:,:,:,:)
      i_nu=me_t*ndim3_nu+ind_nu-1
c ------------ Static part of the vertex --------------------------
      allocate(ktau0(n,n,nsp))
      do ispin=1,nsp
        ktau0(:,:,ispin)=k_pw(:,:,0,1,ispin)
        if(ll==2) ktau0(:,:,ispin)=ktau0(:,:,ispin)+k_pw(:,:,0,2,ispin)
      enddo
      lambda_stat=(0.d0,0.d0)
      call vertex_vk_loc(ktau0,vloc,lambda_stat,nsp,fif,n,np)
      deallocate(ktau0)
c ------------ Dynamic part of the vertex + Q ---------------------
      allocate(lambda_tau(n,n,ll,2,0:n_tau/2,nsp))
      lambda_tau=(0.d0,0.d0)
      if(spin_flips) then
        allocate(qt_loc(np,np,ll,2,0:n_tau/2))
        qt_loc=(0.d0,0.d0)
      endif
      do i0_tau=0,n_tau/2
        do ispin=1,nsp
          isp=ispin
          jsp=ispin
          allocate(kloc(n,n,ll,2))
          do il=1,ll
            do it=1,2
              i_tau=i0_tau
              if(it==2) i_tau=n_tau-i_tau
              kloc(:,:,il,it)=k_pw(:,:,i_tau,il,ispin)
            enddo
          enddo
          call vrt_loop_loc(spin_flips,i0_tau,isp,jsp,ll,kloc,
     &                      w_tau(1,1,i0_tau),fif,
     &                      lambda_tau(1,1,1,1,i0_tau,ispin),
     &                      n,np,gtau(1,1,1,0,ispin))
          deallocate(kloc)
        enddo
      enddo  !! over i0_tau
      if(spin_flips) then
        if(nspin==1.and.irel/=2) qt_loc=2.d0*qt_loc
c ----- T calculation ----------------------------
        call t_loc(qt_loc,vloc,ind_nu,nomnu,wloc,np)
c ----------------- GT part of the vertex -----------------------------
        do i0_tau=0,n_tau/2
          do it=1,2
            call vertex_gt_loc(it,i0_tau,ll,qt_loc(1,1,1,1,i0_tau),
     &                         fif,lambda_tau,n,np,gtau)
          enddo
        enddo
        deallocate(qt_loc)
      endif
c --------------------------------------------------------------------
      do ispin=1,nsp
        isp=ispin
        jsp=ispin
        call vertex_tau_to_omega_loc(lambda_tau(1,1,1,1,0,ispin),
     &                               lambda_dyn(1,1,1,1,ispin),
     &                               ind_nu,nom,isp,jsp,ll,n)
      enddo
      deallocate(lambda_tau)
      end
