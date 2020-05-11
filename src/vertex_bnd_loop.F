      subroutine vertex_bnd_loop(spin_flips,iq,ind_nu,it_vrt,ll,s2,nom,
     &                           nomnu,k_pw,fif,z_red,a_red,phf,ax,v4_r,
     &                           vmi,vim,vii,w4_r,wmi,wim,wii,v2_q,
     &                           lambda_stat,lambda_dyn,nsp,mode,nrr3,
     &                           indrr3)
c     mode = 0 - for regular
c     mode = 1 - for XY
c     mode = 2 - for Sigma_2_14
c     mode = 3 - for Sigma_2_56 with static W
c     mode = 4 - for Sigma_2_56 with dynamic W
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
      integer, intent(in) :: iq,ind_nu,nom,nomnu,s2,it_vrt,ll,nsp,mode,
     &                       nrr3,
     &                       indrr3(nr_full_red,nr_full_red,nqdiv_c)
      real*8, intent(in) :: v4_r(n_pbmt_red,n_pbmt_red,nqdiv_c),
     &                      vmi(n_pbmt_red,nr_full_red,nqdiv_c),
     &                      vim(nr_full_red,n_pbmt_red,nqdiv_c),
     &                      vii(nrr3),
     &                      w4_r(n_pbmt_red,n_pbmt_red,nqdiv_c,
     &                           0:n_tau/2),
     &                      wmi(n_pbmt_red,nr_full_red,nqdiv_c,
     &                          0:n_tau/2),
     &                      wim(nr_full_red,n_pbmt_red,nqdiv_c,
     &                          0:n_tau/2),
     &                      wii(nrr3,0:n_tau/2),
     &                      fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,ll,
     &                               nqdiv_c,nsp),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,
     &                                nspin),
     &                          phf(nr_full_red,nqdiv_c),
     &                          v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,
     &                             nspin)
      complex*16, intent(out) :: lambda_stat(nbndf_bnd,nbndf_bnd,
     &                                       nqdiv_c,nsp),
     &                           lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,
     &                                      nqdiv_c,nsp) 
      integer :: i0_tau,ispin,it,i_nu,isp,jsp,k,il,i_tau
      complex*16 :: cdum
      complex*16, allocatable :: lambda_tau(:,:,:,:,:,:,:),
     &                           ktau0(:,:,:,:),kcom(:,:,:,:,:)
      i_nu=me_t*ndim3_nu+ind_nu-1
      if(mode<3) then
c ------------ Static part of the vertex --------------------------
        allocate(ktau0(nbndf_bnd,nbndf_bnd,nqdiv_c,nsp))
        do ispin=1,nsp
          do k=1,nqdiv_c
            ktau0(:,:,k,ispin)=k_pw(:,:,0,1,k,ispin)
            if(ll==2) ktau0(:,:,k,ispin)=ktau0(:,:,k,ispin)
     &                                  +k_pw(:,:,0,2,k,ispin)
          enddo
        enddo
        lambda_stat=(0.d0,0.d0)
        call vertex_vk(iq,ktau0,z_red,a_red,phf,ax,v4_r,vmi,vim,vii,
     &                 lambda_stat,nsp,mode,fif,0,nrr3,indrr3)
        deallocate(ktau0)
c ------------ Dynamic part of the vertex + Q ---------------------
        allocate(lambda_tau(nbndf_bnd,nbndf_bnd,nqdiv_c,ll,2,0:n_tau/2,
     &                      nsp))
        lambda_tau=(0.d0,0.d0)
      endif
      if(spin_flips) then
        allocate(qt_pw(n_pbtot_red,n_pbtot_red,nqdiv_c,ll,2,0:n_tau/2))
        qt_pw=(0.d0,0.d0)
      endif
      do i0_tau=0,n_tau/2
        do ispin=1,nsp
          if(mod(mode,2)==0) then
            isp=ispin
            jsp=ispin
          else
            isp=nspin
            jsp=1
          endif
          if(mode<3) then
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
            call vrt_loop_mm(spin_flips,i0_tau,isp,jsp,iq,ll,
     &                       kcom,z_red,w4_r(1,1,1,i0_tau),fif,
     &                       lambda_tau(1,1,1,1,1,i0_tau,ispin),mode,
     &                       it_vrt)
            call vrt_loop_im(spin_flips,i0_tau,isp,jsp,iq,ll,
     &                       kcom,z_red,a_red,phf,ax,wim(1,1,1,i0_tau),
     &                       fif,lambda_tau(1,1,1,1,1,i0_tau,ispin),
     &                       mode,it_vrt)
            deallocate(kcom)
c ------------- Common for MI and II ----------------------------------
            allocate(kcom(nbndf_bnd,nrel*nr_full_red,nqdiv_c,ll,2))
            do il=1,ll
              do it=1,2
                i_tau=i0_tau
                if(it==2) i_tau=n_tau-i_tau
                do k=1,nqdiv_c
                  call k_real_space_mi_ii(iq,k,
     &                                    k_pw(1,1,i_tau,il,k,ispin),
     &                                    kcom(1,1,k,il,it),a_red,phf,
     &                                    isp,jsp)
                enddo  !! over k
              enddo
            enddo
            call vrt_loop_mi(spin_flips,i0_tau,isp,jsp,iq,ll,
     &                       kcom,z_red,ax,phf,wmi(1,1,1,i0_tau),fif,
     &                       lambda_tau(1,1,1,1,1,i0_tau,ispin),
     &                       mode,it_vrt)
            call vrt_loop_ii(spin_flips,i0_tau,isp,jsp,iq,ll,
     &                       kcom,a_red,phf,ax,wii(1,i0_tau),
     &                       lambda_tau(1,1,1,1,1,i0_tau,ispin),
     &                       mode,it_vrt,nrr3,indrr3)
            deallocate(kcom)
            call vertex_wk_sing(iq,i0_tau,isp,jsp,ll,
     &                          k_pw(1,1,0,1,1,ispin),
     &                          lambda_tau(1,1,1,1,1,i0_tau,ispin))
          else
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
            call vrt_loop_mm(spin_flips,i0_tau,isp,jsp,iq,ll,
     &                       kcom,z_red,w4_r(1,1,1,i0_tau),fif,cdum,
     &                       mode,it_vrt)
            call vrt_loop_im(spin_flips,i0_tau,isp,jsp,iq,ll,
     &                       kcom,z_red,a_red,phf,ax,wim(1,1,1,i0_tau),
     &                       fif,cdum,mode,it_vrt)
            deallocate(kcom)
c ------------- Common for MI and II ----------------------------------
            allocate(kcom(nbndf_bnd,nrel*nr_full_red,nqdiv_c,ll,2))
            do il=1,ll
              do it=1,2
                i_tau=i0_tau
                if(it==2) i_tau=n_tau-i_tau
                do k=1,nqdiv_c
                  call k_real_space_mi_ii(iq,k,
     &                                    k_pw(1,1,i_tau,il,k,ispin),
     &                                    kcom(1,1,k,il,it),a_red,phf,
     &                                    isp,jsp)
                enddo  !! over k
              enddo
            enddo
            call vrt_loop_mi(spin_flips,i0_tau,isp,jsp,iq,ll,
     &                       kcom,z_red,ax,phf,wmi(1,1,1,i0_tau),fif,
     &                       cdum,mode,it_vrt)
            call vrt_loop_ii(spin_flips,i0_tau,isp,jsp,iq,ll,
     &                       kcom,a_red,phf,ax,wii(1,i0_tau),
     &                       cdum,mode,it_vrt,nrr3,indrr3)
            deallocate(kcom)
          endif
        enddo
      enddo  !! over i0_tau
      if(spin_flips) then
        if(nspin==1.and.irel/=2) qt_pw=2.d0*qt_pw
c ----- T calculation ----------------------------
        if(ll==1) then
          call t_reciprocal_space_static(qt_pw,v2_q,iq)
        else if(ll==2) then
          if(mode==3) then
            call t_recipr_tau(qt_pw,v2_q,iq)
          else
            call t_q_space(qt_pw,v2_q,ind_nu,iq,nomnu,wnu_vrt)
          endif
        endif
        if(iq==1.and.i_nu==0.and.s2==1.and.it_vrt==1)
     &            call timel('**** T_RECIPR finished *************')
c ----------------- GT part of the vertex -----------------------------
        do i0_tau=0,n_tau/2
          do it=1,2
            if(mode<3) then
              call vertex_gt(it,i0_tau,iq,ll,qt_pw(1,1,1,1,1,i0_tau),
     &                       fif,z_red,phf,ax,lambda_tau,mode)
            else
              call vertex_gt(it,i0_tau,iq,ll,qt_pw(1,1,1,1,1,i0_tau),
     &                       fif,z_red,phf,ax,cdum,mode)
            endif
          enddo
        enddo
        deallocate(qt_pw)
        if(iq==1.and.i_nu==0.and.s2==1.and.it_vrt==1)
     &            call timel('**** VERTEX_GT finished ************')
      endif
c --------------------------------------------------------------------
      if(mode>2.or.(mode==2.and.it_vrt==2)) goto 1
      do ispin=1,nsp
        if(mod(mode,2)==0) then
          isp=ispin
          jsp=ispin
        else
          isp=nspin
          jsp=1
        endif
        call vertex_tau_to_omega_bnd_par(lambda_tau(1,1,1,1,1,0,ispin),
     &                                   lambda_dyn(1,1,1,1,1,ispin),
     &                                   ind_nu,iq,nom,isp,jsp,ll)
      enddo
      if(iq==1.and.i_nu==0.and.s2==1.and.it_vrt==1)
     &          call timel('**** VERTEX_TAU_to_OM finished *****')
1     if(mode<3) deallocate(lambda_tau)
      end
