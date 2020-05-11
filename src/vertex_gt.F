      subroutine vertex_gt(it,i0_tau,iq,ll,q_pw,fif,z_red,phf,ax,
     &                     lambda_tau,mode)
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
      integer, intent(in) :: it,i0_tau,iq,ll,mode
      real*8, intent(in) :: fif(nrel*nindm_fif,maxel_red,maxel_red,
     &                          nsort,nspin)
      complex*16, intent(in) :: q_pw(n_pbtot_red,n_pbtot_red,nqdiv_c,ll,
     &                               2),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          phf(nr_full_red,nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,
     &                             nspin)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        nqdiv_c,ll,2,0:n_tau/2,
     &                                        nspin)
      integer :: i_tau,ispin,il,k
      complex*16, allocatable :: vx(:,:,:,:,:),t_rs(:,:,:,:)
      i_tau=i0_tau
      if(it==2) i_tau=n_tau-i0_tau
c -------------------------------- MM --------------------------------
      allocate(t_rs(n_pbmt_red,n_pbmt_red,nqdiv_c,ll))
      call t_rs_mm(it,ll,q_pw,t_rs,n_pbtot_red)
      allocate(vx(nfun_red,nfun_red,nqdiv_c,ll,nspin))
      vx=(0.d0,0.d0)
      if(irel/=2) then
        do ispin=1,nspin
          call vertex_mm_gt_r(i_tau,ll,t_rs,g_rs_mm(1,1,1,0,1,ispin),
     &                        vx(1,1,1,1,ispin),fif(1,1,1,1,ispin))
        enddo
      else
        call vertex_mm_gt_c(i_tau,ll,t_rs,g_rs_mm(1,1,1,0,1,1),vx,
     &                      fif(1,1,1,1,1))
      endif
      deallocate(t_rs)
      do ispin=1,nspin
        do il=1,ll
          call vertex_rr_to_k_all(vx(1,1,1,il,ispin),nfun_red**2)
          if(mode>2) then
            do k=1,nqdiv_c
              vdy_store(1:nfun_red,1:nfun_red,i_tau,il,k,ispin)=
     &                      vx(:,:,k,il,ispin)
            enddo
          endif
          if(mode<3) then
            call vertex_rs_to_bnd_mm(vx(1,1,1,il,ispin),
     &                             lambda_tau(1,1,1,il,it,i0_tau,ispin),
     &                               iq,z_red,ispin,ispin)
          endif
        enddo
      enddo
      deallocate(vx)
c -------------------------------- MI --------------------------------
      allocate(t_rs(nr_full_red,n_pbmt_red,nqdiv_c,ll)) 
      call t_rs_im(it,ll,q_pw,t_rs,phf)
      allocate(vx(nfun_red,nrel*nr_full_red,nqdiv_c,ll,nspin))
      vx=(0.d0,0.d0)
      if(irel/=2) then
        do ispin=1,nspin 
          call vertex_mi_gt_r(i_tau,ll,t_rs,g_rs_mi(1,1,1,0,1,ispin),
     &                        vx(1,1,1,1,ispin),fif(1,1,1,1,ispin))
        enddo
      else 
        call vertex_mi_gt_c(i_tau,ll,t_rs,g_rs_mi(1,1,1,0,1,1),vx,
     &                      fif(1,1,1,1,1))
      endif
      deallocate(t_rs)
      do ispin=1,nspin
        do il=1,ll
          call vertex_rr_to_k_all(vx(1,1,1,il,ispin),
     &                            nfun_red*nrel*nr_full_red)
          if(mode>2) then
            do k=1,nqdiv_c
              vdy_store(1:nfun_red,nfun_red+1:nrs_red,i_tau,il,k,ispin)=
     &                    vx(:,:,k,il,ispin)
            enddo
          endif
          if(mode<3) then
            call vertex_rs_to_bnd_mi(vx(1,1,1,il,ispin),
     &                             lambda_tau(1,1,1,il,it,i0_tau,ispin),
     &                               iq,z_red,ax,phf,ispin,ispin)
          endif
        enddo
      enddo
      deallocate(vx)
c -------------------------------- IM --------------------------------
      allocate(t_rs(n_pbmt_red,nr_full_red,ll,nqdiv_c)) 
      call t_rs_mi(it,ll,q_pw,t_rs,iq,phf)
      allocate(vx(nrel*nr_full_red,nfun_red,nqdiv_c,ll,nspin))
      vx=(0.d0,0.d0)
      if(irel/=2) then
        do ispin=1,nspin 
          call vertex_im_gt_r(i_tau,ll,t_rs,g_rs_im(1,1,1,0,1,ispin),
     &                          vx(1,1,1,1,ispin),fif(1,1,1,1,ispin))
        enddo
      else 
        call vertex_im_gt_c(i_tau,ll,t_rs,g_rs_im(1,1,1,0,1,1),vx,fif)
      endif
      deallocate(t_rs)
      do ispin=1,nspin
        do il=1,ll
          call vertex_rr_to_k_all(vx(1,1,1,il,ispin),
     &                            nfun_red*nrel*nr_full_red)
          if(mode>2) then
            do k=1,nqdiv_c
              vdy_store(nfun_red+1:nrs_red,1:nfun_red,i_tau,il,k,ispin)=
     &                    vx(:,:,k,il,ispin)
            enddo
          endif
          if(mode<3) then
            call vertex_rs_to_bnd_im(vx(1,1,1,il,ispin),
     &                             lambda_tau(1,1,1,il,it,i0_tau,ispin),
     &                               iq,z_red,ax,phf,ispin,ispin)
          endif
       enddo
      enddo
      deallocate(vx)
c -------------------------------- II --------------------------------
      allocate(t_rs(nr_full_red,nr_full_red,ll,nqdiv_c))  
      call t_rs_ii(it,ll,q_pw,t_rs,phf,iq)
      allocate(vx(nrel*nr_full_red,nrel*nr_full_red,nqdiv_c,ll,nspin))
      vx=(0.d0,0.d0)
      if(irel/=2) then
        do ispin=1,nspin 
          call vertex_ii_gt_r(i_tau,ll,t_rs,g_rs_ii(1,1,1,0,1,ispin),
     &                        vx(1,1,1,1,ispin))
        enddo
      else 
        call vertex_ii_gt_c(i_tau,ll,t_rs,g_rs_ii(1,1,1,0,1,1),vx)
      endif
      deallocate(t_rs)
      do ispin=1,nspin
        do il=1,ll
          call vertex_rr_to_k_all(vx(1,1,1,il,ispin),
     &                            (nrel*nr_full_red)**2)
          if(mode>2) then
            do k=1,nqdiv_c
              vdy_store(nfun_red+1:nrs_red,nfun_red+1:nrs_red,i_tau,
     &                  il,k,ispin)=vx(:,:,k,il,ispin)
            enddo
          endif
          if(mode<3) then
            call vertex_rs_to_bnd_ii(vx(1,1,1,il,ispin),
     &                             lambda_tau(1,1,1,il,it,i0_tau,ispin),
     &                               iq,ax,phf,ispin,ispin)
          endif
        enddo
      enddo
      deallocate(vx)
      end
