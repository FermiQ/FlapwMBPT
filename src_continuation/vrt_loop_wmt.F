      subroutine vrt_loop_wmt(spin_flips,i0_tau,isp,jsp,iq,ll,kcom,
     &                        z_red,w_tau,fif,lambda_tau)
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
      integer, intent(in) :: i0_tau,iq,ll,isp,jsp
      real*8, intent(in) :: w_tau(n_pbmtm_red,n_pbmtm_red,natom),
     &                      fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: kcom(nbndf_bnd,nfun_red,nqdiv_c,ll,2),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        nqdiv_c,ll,2)
      integer :: il,it,i_tau,k,iatom,isort,n,i0
      complex*16, allocatable :: k_rs(:,:,:,:,:),vx(:,:,:),
     &                           kloc(:,:,:,:,:)
      allocate(kloc(maxel_red,maxel_red,natom,ll,2))
      if(spin_flips) then
        allocate(k_rs(nfun_red,nfun_red,nqdiv_c,ll,2))
        do il=1,ll
          do it=1,2
            do k=1,nqdiv_c
              call k_real_space_mm(k,kcom(1,1,k,il,it),
     &                             k_rs(1,1,k,il,it),z_red,isp)
            enddo  !! over k
            call k_k_to_rr_all(k_rs(1,1,1,il,it),nfun_red**2)
            do iatom=1,natom
              i0=io_lem_red(iatom)
              isort=is(iatom)
              n=lfunm_red(isort)
              kloc(1:n,1:n,iatom,il,it)=
     &            k_rs(i0:i0+n-1,i0:i0+n-1,1,il,it)
            enddo
          enddo
        enddo
      else
        do il=1,ll
          do it=1,2
            call k_real_space_loc(kcom(1,1,1,il,it),kloc(1,1,1,il,it),
     &                            z_red,isp)
          enddo
        enddo
      endif
      allocate(vx(maxel_red,maxel_red,natom))
      do il=1,ll
        do it=1,2
          vx=(0.d0,0.d0)
          if(irel<2) call vertex_wk_loc_r(w_tau,kloc(1,1,1,il,it),vx,
     &                                    fif,isp,jsp)
          if(irel==2) call vertex_wk_loc_c(w_tau,kloc(1,1,1,il,it),vx,
     &                                     fif)
          call vertex_loc_to_bnd(vx,lambda_tau(1,1,1,il,it),iq,z_red,
     &                           isp,jsp)
          if(spin_flips) then
            i_tau=i0_tau
            if(it==2) i_tau=n_tau-i_tau
            if(irel<2) call q_tau_mm_r(i_tau,ll,it,il,
     &                                 qt_pw(1,1,1,il,it,i0_tau),
     &                                 g_rs_mm(1,1,1,0,1,isp),k_rs,
     &                                 fif(1,1,1,1,isp),n_pbmt_red)
            if(irel==2) call q_tau_mm_c(i_tau,ll,it,il,
     &                                  qt_pw(1,1,1,il,it,i0_tau),
     &                                  g_rs_mm(1,1,1,0,1,isp),k_rs,
     &                                  fif(1,1,1,1,isp),n_pbmt_red)
          endif
        enddo
      enddo
      deallocate(kloc)
      if(spin_flips) deallocate(k_rs)
      deallocate(vx)
      end