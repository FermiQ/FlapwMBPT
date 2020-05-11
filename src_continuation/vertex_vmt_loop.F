      subroutine vertex_vmt_loop(spin_flips,iq,ind_nu,it_vrt,ll,s2,nom,
     &                           nomnu,k_pw,fif,z_red,vloc,
     &                           lambda_stat,lambda_dyn,nsp)
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
      integer, intent(in) :: iq,ind_nu,nom,nomnu,s2,it_vrt,ll,nsp
      real*8, intent(in) :: vloc(n_pbmtm_red,n_pbmtm_red,natom),
     &                      fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,ll,
     &                               nqdiv_c,nsp),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(out) :: lambda_stat(nbndf_bnd,nbndf_bnd,
     &                                       nqdiv_c,nsp),
     &                           lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,
     &                                      nqdiv_c,nsp)
      integer :: i0_tau,ispin,it,i_nu,isp,jsp,k,il,i_tau
      complex*16, allocatable :: lambda_tau(:,:,:,:,:,:,:),
     &                           ktau0(:,:,:,:),kcom(:,:,:,:,:)
      i_nu=me_t*ndim3_nu+ind_nu-1
c ------------ Static part of the vertex --------------------------
      allocate(ktau0(nbndf_bnd,nbndf_bnd,nqdiv_c,nsp))
      do ispin=1,nsp
        do k=1,nqdiv_c
          ktau0(:,:,k,ispin)=k_pw(:,:,0,1,k,ispin)
          if(ll==2) ktau0(:,:,k,ispin)=ktau0(:,:,k,ispin)
     &                                +k_pw(:,:,0,2,k,ispin)
        enddo
      enddo
      lambda_stat=(0.d0,0.d0)
      call vertex_vk_wmt(iq,ktau0,z_red,vloc,lambda_stat,nsp,fif)
      deallocate(ktau0)
c ------------ Dynamic part of the vertex (from Q only) ------------
      if(spin_flips) then
        allocate(lambda_tau(nbndf_bnd,nbndf_bnd,nqdiv_c,ll,2,
     &                      0:n_tau/2,nsp))
        lambda_tau=(0.d0,0.d0)
        allocate(qt_pw(n_pbmt_red,n_pbmt_red,nqdiv_c,ll,2,0:n_tau/2))
        qt_pw=(0.d0,0.d0)
        do i0_tau=0,n_tau/2
          do ispin=1,nsp
            isp=ispin
            jsp=ispin
c ------------- Common for MM and IM ----------------------------------
            allocate(kcom(nbndf_bnd,nfun_red,nqdiv_c,ll,2))
            do il=1,ll
              do it=1,2
                i_tau=i0_tau
                if(it==2) i_tau=n_tau-i_tau
                do k=1,nqdiv_c
                  call k_real_space_mm_im(iq,k,
     &                                    k_pw(1,1,i_tau,il,k,ispin),
     &                                    kcom(1,1,k,il,it),z_red,isp,
     &                                    jsp)
                enddo  !! over k
              enddo
            enddo
            call vrt_loop_vmt(spin_flips,i0_tau,isp,jsp,iq,ll,kcom,
     &                        z_red,fif)
            deallocate(kcom)
          enddo
        enddo  !! over i0_tau
        if(nspin==1.and.irel/=2) qt_pw=2.d0*qt_pw
c ----- T calculation ----------------------------
        call t_rr_stat(qt_pw,vloc,ind_nu,nomnu)
        if(iq==1.and.i_nu==0.and.s2==1.and.it_vrt==1)
     &            call timel('**** T_RECIPR finished *************')
c ----------------- GT part of the vertex -----------------------------
        do i0_tau=0,n_tau/2
          do it=1,2
            call vertex_gt_wmt(it,i0_tau,iq,ll,qt_pw(1,1,1,1,1,i0_tau),
     &                         fif,z_red,lambda_tau)
          enddo
        enddo
        deallocate(qt_pw)
        if(iq==1.and.i_nu==0.and.s2==1.and.it_vrt==1)
     &            call timel('**** VERTEX_GT finished ************')
c --------------------------------------------------------------------
        do ispin=1,nsp
          isp=ispin
          jsp=ispin
         call vertex_tau_to_omega_bnd_par(lambda_tau(1,1,1,1,1,0,ispin),
     &                                    lambda_dyn(1,1,1,1,1,ispin),
     &                                    ind_nu,iq,nom,isp,jsp,ll)
        enddo
        if(iq==1.and.i_nu==0.and.s2==1.and.it_vrt==1)
     &          call timel('**** VERTEX_TAU_to_OM finished *****')
        deallocate(lambda_tau)
      endif
      end
