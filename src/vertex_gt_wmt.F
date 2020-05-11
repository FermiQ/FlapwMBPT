      subroutine vertex_gt_wmt(it,i0_tau,iq,ll,q_pw,fif,z_red,
     &                         lambda_tau)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,i0_tau,iq,ll
      real*8, intent(in) :: fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: q_pw(n_pbmt_red,n_pbmt_red,nqdiv_c,ll,
     &                               2),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        nqdiv_c,ll,2,0:n_tau/2,
     &                                        nspin)
      integer :: i_tau,ispin,il
      complex*16, allocatable :: vx(:,:,:,:,:)
      i_tau=i0_tau
      if(it==2) i_tau=n_tau-i0_tau
c -------------------------------- MM --------------------------------
      allocate(vx(nfun_red,nfun_red,nqdiv_c,ll,nspin))
      vx=(0.d0,0.d0)
      if(irel/=2) then
        do ispin=1,nspin
          call vertex_mm_gt_r(i_tau,ll,q_pw(1,1,1,1,it),
     &                        g_rs_mm(1,1,1,0,1,ispin),
     &                        vx(1,1,1,1,ispin),fif(1,1,1,1,ispin))
        enddo
      else
        call vertex_mm_gt_c(i_tau,ll,q_pw(1,1,1,1,it),
     &                      g_rs_mm(1,1,1,0,1,1),vx,fif(1,1,1,1,1))
      endif
      do ispin=1,nspin
        do il=1,ll
          call vertex_rr_to_k_all(vx(1,1,1,il,ispin),nfun_red**2)
          call vertex_rs_to_bnd_mm(vx(1,1,1,il,ispin),
     &                             lambda_tau(1,1,1,il,it,i0_tau,ispin),
     &                             iq,z_red,ispin,ispin)
        enddo
      enddo
      deallocate(vx)
      end
