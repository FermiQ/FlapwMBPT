      subroutine vrt_loop_vmt(spin_flips,i0_tau,isp,jsp,iq,ll,kcom,
     &                        z_red,fif)
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
      real*8, intent(in) :: fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: kcom(nbndf_bnd,nfun_red,nqdiv_c,ll,2),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      integer :: il,it,i_tau,k
      complex*16, allocatable :: k_rs(:,:,:,:,:)
      allocate(k_rs(nfun_red,nfun_red,nqdiv_c,ll,2))
      do il=1,ll
        do it=1,2
          do k=1,nqdiv_c
            call k_real_space_mm(k,kcom(1,1,k,il,it),
     &                           k_rs(1,1,k,il,it),z_red,isp)
          enddo  !! over k
          call k_k_to_rr_all(k_rs(1,1,1,il,it),nfun_red**2)
        enddo
      enddo
      do il=1,ll
        do it=1,2
          i_tau=i0_tau
          if(it==2) i_tau=n_tau-i_tau
          if(irel<2) call q_tau_mm_r(i_tau,ll,it,il,
     &                               qt_pw(1,1,1,il,it,i0_tau),
     &                               g_rs_mm(1,1,1,0,1,isp),k_rs,
     &                               fif(1,1,1,1,isp),n_pbmt_red)
          if(irel==2) call q_tau_mm_c(i_tau,ll,it,il,
     &                                qt_pw(1,1,1,il,it,i0_tau),
     &                                g_rs_mm(1,1,1,0,1,isp),k_rs,
     &                                fif(1,1,1,1,isp),n_pbmt_red)
        enddo
      enddo
      deallocate(k_rs)
      end
