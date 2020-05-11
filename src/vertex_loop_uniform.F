      subroutine vertex_loop_uniform(spin_flips,it_vrt,k_pw,fif,
     &                               z_red,a_red,phf,ax,v4_r,vmi,vim,
     &                               vii,w4_r,wmi,wim,wii,v2_q,
     &                               lambda_stat,lambda_dyn,nsp,nrr3,
     &                               indrr3)
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
      integer, intent(in) :: it_vrt,nsp,nrr3,
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
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,
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
     &                           lambda_dyn(nbndf_bnd,nbndf_bnd,
     &                                      0:n_omega,nqdiv_c,nsp) 
      integer :: i0_tau,ispin,it,isp,jsp,k,i_tau,ind_tau
      complex*16, allocatable :: lambda_tau(:,:,:,:,:,:),
     &                           ktau0(:,:,:,:),kcom(:,:,:,:)
c ------------ Static part of the vertex --------------------------
      allocate(ktau0(nbndf_bnd,nbndf_bnd,nqdiv_c,nsp))
      do ispin=1,nsp
        do k=1,nqdiv_c
          ktau0(:,:,k,ispin)=k_pw(:,:,0,k,ispin)
        enddo
      enddo
      lambda_stat=(0.d0,0.d0)
      call vertex_vk(1,ktau0,z_red,a_red,phf,ax,v4_r,vmi,vim,vii,
     &               lambda_stat,nsp,0,fif,0,nrr3,indrr3)
      deallocate(ktau0)
c ------------ Dynamic part of the vertex + Q ---------------------
      allocate(lambda_tau(nbndf_bnd,nbndf_bnd,nqdiv_c,2,0:n_tau/2,nsp))
      lambda_tau=(0.d0,0.d0)
      if(spin_flips) then
        allocate(qt_pw(n_pbtot_red,n_pbtot_red,nqdiv_c,1,2,0:n_tau/2))
        qt_pw=(0.d0,0.d0)
      endif
      do ind_tau=1,ndim3_tau
        i0_tau=me_t*ndim3_tau+ind_tau-1
        do ispin=1,nsp
          isp=ispin
          jsp=ispin
c ------------- Common for MM and IM -------------------------------
          allocate(kcom(nbndf_bnd,nfun_red,nqdiv_c,2))
          do it=1,2
            i_tau=i0_tau
            if(it==2) i_tau=n_tau-i_tau
            do k=1,nqdiv_c
              call k_real_space_mm_im(1,k,k_pw(1,1,i_tau,k,ispin),
     &                                kcom(1,1,k,it),z_red,isp,jsp)
            enddo  !! over k
          enddo
          call vrt_loop_mm(spin_flips,i0_tau,isp,jsp,1,1,kcom,z_red,
     &                     w4_r(1,1,1,i0_tau),fif,
     &                     lambda_tau(1,1,1,1,i0_tau,ispin),0,it_vrt)
          call vrt_loop_im(spin_flips,i0_tau,isp,jsp,1,1,kcom,z_red,
     &                     a_red,phf,ax,wim(1,1,1,i0_tau),fif,
     &                     lambda_tau(1,1,1,1,i0_tau,ispin),0,it_vrt)
          deallocate(kcom)
c ------------- Common for MI and II ----------------------------------
          allocate(kcom(nbndf_bnd,nrel*nr_full_red,nqdiv_c,2))
          do it=1,2
            i_tau=i0_tau
            if(it==2) i_tau=n_tau-i_tau
            do k=1,nqdiv_c
              call k_real_space_mi_ii(1,k,k_pw(1,1,i_tau,k,ispin),
     &                                kcom(1,1,k,it),a_red,phf,isp,jsp)
            enddo  !! over k
          enddo
          call vrt_loop_mi(spin_flips,i0_tau,isp,jsp,1,1,kcom,z_red,
     &                     ax,phf,wmi(1,1,1,i0_tau),fif,
     &                     lambda_tau(1,1,1,1,i0_tau,ispin),0,it_vrt)
          call vrt_loop_ii(spin_flips,i0_tau,isp,jsp,1,1,kcom,a_red,phf,
     &                     ax,wii(1,i0_tau),
     &                     lambda_tau(1,1,1,1,i0_tau,ispin),0,it_vrt,
     &                     nrr3,indrr3)
          deallocate(kcom)
        enddo
      enddo  !! over ind_tau
      if(spin_flips) then
        if(nspin==1.and.irel/=2) qt_pw=2.d0*qt_pw
        if(nproc_t/=1) call DGOP(qt_pw,
     &                             2*n_pbtot_red**2*nqdiv_c*(n_tau+2),
     &                             '  +',comm_t)
c ----- T calculation ----------------------------
        call t_reciprocal_space_static(qt_pw,v2_q,1)
        if(it_vrt==1) call timel('**** T_RECIPR finished *************')
c ----------------- GT part of the vertex -----------------------------
        do ind_tau=1,ndim3_tau
          i0_tau=me_t*ndim3_tau+ind_tau-1
          do it=1,2
            call vertex_gt(it,i0_tau,1,1,qt_pw(1,1,1,1,1,i0_tau),fif,
     &                     z_red,phf,ax,lambda_tau,0)
          enddo
        enddo
        deallocate(qt_pw)
        if(it_vrt==1)call timel('**** VERTEX_GT finished ************')
      endif
	  if(nproc_t/=1) call DGOP(lambda_tau,
     &                           2*nbndf_bnd**2*nqdiv_c*(n_tau+2)*nsp,
     &                           '  +',comm_t)
c --------------------------------------------------------------------
      do ispin=1,nsp
        isp=ispin
        jsp=ispin
        call vertex_tau_to_omega_uniform(lambda_tau(1,1,1,1,0,ispin),
     &                                   lambda_dyn(1,1,0,1,ispin),
     &                                   isp,jsp)
      enddo
      if(it_vrt==1)call timel('**** VERTEX_TAU_to_OM finished *****')
      deallocate(lambda_tau)
      end
