      subroutine vertex_vk_wmt(iq,k_pw,z_red,vloc,lambda_stat,nsp,fif)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq,nsp
      real*8, intent(in) :: vloc(n_pbmtm_red,n_pbmtm_red,natom),
     &                      fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,nqdiv_c,nsp),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(out) :: lambda_stat(nbndf_bnd,nbndf_bnd,
     &                                       nqdiv_c,nsp)
      integer :: ispin,k,isp,jsp,nrr
      complex*16, allocatable :: krs0(:,:,:),vx(:,:,:),kcom(:,:,:)
      nrr=nrel*nr_full_red
      do ispin=1,nsp
        isp=ispin
        jsp=ispin
c ------------- Common for MM and IM ----------------------------------
        allocate(kcom(nbndf_bnd,nfun_red,nqdiv_c))
        do k=1,nqdiv_c
          call k_real_space_mm_im(iq,k,k_pw(1,1,k,ispin),kcom(1,1,k),
     &                            z_red,isp,jsp)
        enddo  !! over k
        allocate(krs0(maxel_red,maxel_red,natom))
        call k_real_space_loc(kcom,krs0,z_red,isp)
        deallocate(kcom)
	    allocate(vx(maxel_red,maxel_red,natom))
        vx=(0.d0,0.d0)
        if(irel<2) call vertex_wk_loc_r(vloc,krs0,vx,fif,isp,jsp)
        if(irel==2) call vertex_wk_loc_c(vloc,krs0,vx,fif)
        deallocate(krs0)
        call vertex_loc_to_bnd(vx,lambda_stat(1,1,1,ispin),iq,z_red,
     &                         isp,jsp)
        deallocate(vx)
      enddo  !! over ispin
      end
