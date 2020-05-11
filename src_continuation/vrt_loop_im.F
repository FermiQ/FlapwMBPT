      subroutine vrt_loop_im(spin_flips,i0_tau,isp,jsp,iq,ll,kcom,z_red,
     &                       a_red,phf,ax,wim,fif,lambda_tau,mode,
     &                       it_vrt)
c     mode = 0 - for regular
c     mode = 1 - for XY
c     mode = 2 - for Sigma_2_14
c     mode = 3 - for or Sigma_2_56 with static W
c     mode = 4 - for or Sigma_2_56 with dynamic W
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
      integer, intent(in) :: i0_tau,iq,ll,isp,jsp,mode,it_vrt
      real*8, intent(in) :: wim(nr_full_red,n_pbmt_red,nqdiv_c),
     &                      fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: kcom(nbndf_bnd,nfun_red,nqdiv_c,ll,2),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,
     &                                nspin),
     &                          phf(nr_full_red,nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,
     &                             nspin)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        nqdiv_c,ll,2)
      integer :: il,it,i_tau,nrr,k
      complex*16, allocatable :: k_rs(:,:,:,:,:),vx(:,:,:)
      nrr=nrel*nr_full_red
      allocate(k_rs(nrr,nfun_red,nqdiv_c,ll,2))
      do il=1,ll
        do it=1,2
          do k=1,nqdiv_c
            call k_real_space_im(iq,k,kcom(1,1,k,il,it),
     &                           k_rs(1,1,k,il,it),a_red,phf,isp)
          enddo  !! over k
          call k_k_to_rr_all(k_rs(1,1,1,il,it),nfun_red*nrr)
        enddo
      enddo
      if(mode<3) allocate(vx(nrr,nfun_red,nqdiv_c))
      do il=1,ll
        do it=1,2
          if(mode<3) then
            vx=(0.d0,0.d0)
            if(irel<2) call vertex_wk_im_r(wim,k_rs(1,1,1,il,it),vx,
     &                                     fif(1,1,1,1,jsp))
            if(irel==2) call vertex_wk_im_c(wim,k_rs(1,1,1,il,it),vx,
     &                                      fif(1,1,1,1,jsp))
            call vertex_rr_to_k_all(vx,nfun_red*nrr)
            if(mode==2) then
              i_tau=i0_tau
              if(it==2) i_tau=n_tau-i_tau
              do k=1,nqdiv_c
                vdy_store(nfun_red+1:nrs_red,1:nfun_red,i_tau,il,k,isp)=
     &                    vx(:,:,k)
              enddo
              if(it_vrt==2) goto 1
            endif
            call vertex_rs_to_bnd_im(vx,lambda_tau(1,1,1,il,it),iq,
     &                               z_red,ax,phf,isp,jsp)
          endif
1         if(spin_flips) then
            i_tau=i0_tau
            if(it==2) i_tau=n_tau-i_tau
            if(irel<2) call q_tau_im_r(i_tau,ll,it,il,iq,
     &                                 qt_pw(1,1,1,il,it,i0_tau),
     &                                 g_rs_mi(1,1,1,0,1,isp),
     &                                 g_rs_im(1,1,1,0,1,isp),k_rs,
     &                                 fif(1,1,1,1,isp),phf)
            if(irel==2) call q_tau_im_c(i_tau,ll,it,il,iq,
     &                                  qt_pw(1,1,1,il,it,i0_tau),
     &                                  g_rs_mi(1,1,1,0,1,isp),
     &                                  g_rs_im(1,1,1,0,1,isp),k_rs,
     &                                  fif(1,1,1,1,isp),phf)
          endif
        enddo
      enddo
      deallocate(k_rs)
      if(mode<3) deallocate(vx)
      end
