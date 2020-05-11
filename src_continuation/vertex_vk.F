      subroutine vertex_vk(iq,k_pw,z_red,a_red,phf,ax,v4_r,vmi,vim,vii,
     &                     lambda_stat,nsp,mode,fif,khbr,nrr3,indrr3)
c     mode = 0 - for regular
c     mode = 1 - for XY
c     mode = 2 - for Sigma_2_14
c     mode = 3 - for or Sigma_2_56 with static W
c     mode = 4 - for or Sigma_2_56 with dynamic W
c     khbr = 0 - singular term V(q=0) is treated separately
c     khbr = 1 - singular term V(q=0) is treated as the rest (Yukava)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq,nsp,mode,khbr,nrr3,
     &                       indrr3(nr_full_red,nr_full_red,nqdiv_c)
      real*8, intent(in) :: v4_r(n_pbmt_red,n_pbmt_red,nqdiv_c),
     &                      vmi(n_pbmt_red,nr_full_red,nqdiv_c),
     &                      vim(nr_full_red,n_pbmt_red,nqdiv_c),
     &                      vii(nrr3),
     &                      fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,nqdiv_c,nsp),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,
     &                                nspin),
     &                          phf(nr_full_red,nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,
     &                             nspin)
      complex*16, intent(out) :: lambda_stat(nbndf_bnd,nbndf_bnd,
     &                                       nqdiv_c,nsp)
      integer :: ispin,k,isp,jsp,nrr
      real*8 :: w
      complex*16, allocatable :: krs0(:,:,:),vx(:,:,:),kcom(:,:,:)
      nrr=nrel*nr_full_red
      do ispin=1,nsp
        if(mod(mode,2)==0) then
          isp=ispin
          jsp=ispin
        else
          isp=1
          jsp=nspin
        endif
c ------------- Common for MM and IM ----------------------------------
        allocate(kcom(nbndf_bnd,nfun_red,nqdiv_c))
        do k=1,nqdiv_c
          call k_real_space_mm_im(iq,k,k_pw(1,1,k,ispin),kcom(1,1,k),
     &                            z_red,isp,jsp)
        enddo  !! over k
c -----------------  MM -----------------------------------------------        
        allocate(krs0(nfun_red,nfun_red,nqdiv_c))
        do k=1,nqdiv_c
          call k_real_space_mm(k,kcom(1,1,k),krs0(1,1,k),z_red,isp)
        enddo  !! over k
        call k_k_to_rr_all(krs0,nfun_red**2)
	    allocate(vx(nfun_red,nfun_red,nqdiv_c))
        vx=(0.d0,0.d0)
        if(irel<2) call vertex_wk_mm_r(v4_r,krs0,vx,fif,isp,jsp)
        if(irel==2) call vertex_wk_mm_c(v4_r,krs0,vx,fif)
        deallocate(krs0)
        call vertex_rr_to_k_all(vx,nfun_red**2)
        if(mode==2) then
          vst_store(1:nfun_red,1:nfun_red,:,ind_nu_st,ind_pbr_st,ispin)=
     &                   vx
        endif
        call vertex_rs_to_bnd_mm(vx,lambda_stat(1,1,1,ispin),iq,z_red,
     &                           isp,jsp)
        deallocate(vx)
c -----------------  IM -----------------------------------------------        
        allocate(krs0(nrr,nfun_red,nqdiv_c))
        do k=1,nqdiv_c
          call k_real_space_im(iq,k,kcom(1,1,k),krs0(1,1,k),a_red,phf,
     &                         isp)
        enddo  !! over k
        call k_k_to_rr_all(krs0,nfun_red*nrr)
	    allocate(vx(nrr,nfun_red,nqdiv_c))
        vx=(0.d0,0.d0)
        if(irel<2) call vertex_wk_im_r(vim,krs0,vx,fif(1,1,1,1,jsp))
        if(irel==2) call vertex_wk_im_c(vim,krs0,vx,fif(1,1,1,1,jsp))
        deallocate(krs0)
        call vertex_rr_to_k_all(vx,nfun_red*nrr)
        if(mode==2) then
          vst_store(nfun_red+1:nrs_red,1:nfun_red,:,ind_nu_st,
     &              ind_pbr_st,ispin)=vx
        endif
        call vertex_rs_to_bnd_im(vx,lambda_stat(1,1,1,ispin),iq,z_red,
     &                           ax,phf,isp,jsp)
        deallocate(vx,kcom)
c ------------- Common for MI and II ----------------------------------
        allocate(kcom(nbndf_bnd,nrel*nr_full_red,nqdiv_c))
        do k=1,nqdiv_c
          call k_real_space_mi_ii(iq,k,k_pw(1,1,k,ispin),kcom(1,1,k),
     &                            a_red,phf,isp,jsp)
        enddo  !! over k
c -----------------  MI -----------------------------------------------        
        allocate(krs0(nfun_red,nrr,nqdiv_c))
        do k=1,nqdiv_c
          call k_real_space_mi(iq,k,kcom(1,1,k),krs0(1,1,k),z_red,isp)
        enddo  !! over k
        call k_k_to_rr_all(krs0,nfun_red*nrr)
	    allocate(vx(nfun_red,nrr,nqdiv_c))
        vx=(0.d0,0.d0)
        if(irel<2) call vertex_wk_mi_r(vmi,krs0,vx,fif,isp)
        if(irel==2) call vertex_wk_mi_c(vmi,krs0,vx,fif)
        deallocate(krs0)
        call vertex_rr_to_k_all(vx,nfun_red*nrr)
        if(mode==2) then
          vst_store(1:nfun_red,nfun_red+1:nrs_red,:,ind_nu_st,
     &              ind_pbr_st,ispin)=vx
        endif
        call vertex_rs_to_bnd_mi(vx,lambda_stat(1,1,1,ispin),iq,z_red,
     &                           ax,phf,isp,jsp)
        deallocate(vx)
c -----------------  II -----------------------------------------------      
        allocate(vx(nrr,nrr,nqdiv_c))
        do k=1,nqdiv_c
          call k_real_space_ii(iq,k,kcom(1,1,k),vx(1,1,k),a_red,phf,isp)
        enddo  !! over k
        call k_k_to_rr_all(vx,nrr*nrr)
        call vertex_wk_ii(vii,vx,nrr3,indrr3)
        call vertex_rr_to_k_all(vx,nrr*nrr)
        if(mode==2) then
          vst_store(nfun_red+1:nrs_red,nfun_red+1:nrs_red,:,ind_nu_st,
     &              ind_pbr_st,ispin)=vx
        endif
        call vertex_rs_to_bnd_ii(vx,lambda_stat(1,1,1,ispin),iq,ax,
     &                           phf,isp,jsp)
        deallocate(vx,kcom)
        if(khbr==0) then
          w=8.d0*pi*q2aver_c/amega/nqdiv_c
          if(key_hybrid==1) w=w*alpha_hbr
          do k=1,nqdiv_c
            lambda_stat(:,:,k,ispin)=lambda_stat(:,:,k,ispin)
     &                              +w*k_pw(:,:,k,ispin)
          enddo  !! over k
        endif
      enddo  !! over ispin
      end
