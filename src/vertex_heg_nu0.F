      subroutine vertex_heg_nu0(iter_vrt,iq,lambda_dyn,lambda_stat,
     &                          conv_table,conv_1,spin_flips,ntr,q)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      logical, intent(in) :: spin_flips
      integer, intent(in) :: iter_vrt,iq,ntr
      real*8, intent(in) :: q(*)
      real*8, intent(inout) :: conv_table(2:ntr,3)
      real*8, intent(out) :: conv_1
      complex*16, intent(out) :: lambda_stat(nrr_div,nspin),
     &                           lambda_dyn(nrr_div,0:n_omega,nspin)
      logical :: conv
      integer :: it,ispin
      real*8 :: conv_0,lam_max
      complex*16, allocatable :: lambda_dyn_tau(:,:,:,:),
     &                           k0_pw(:,:,:,:),k_pw(:,:,:,:),
     &                           q_tau(:,:),t_tau(:,:),
     &                           vrt_dyn_old(:,:,:),
     &                           vrt_stat_old(:,:)
      conv_1=0.d0
c --------------------------------------------------------------------
      lambda_stat=(0.d0,0.d0)
      lambda_dyn=(0.d0,0.d0)
	allocate(k0_pw(nrr_div,2,ndim3_tau,nspin))
	allocate(k_pw(nrr_div,2,ndim3_tau,nspin))
	allocate(lambda_dyn_tau(nrr_div,2,ndim3_tau,nspin))
	allocate(vrt_dyn_old(nrr_div,0:n_omega,nspin))
	allocate(vrt_stat_old(nrr_div,nspin))
      if(spin_flips.and.iter_vrt>1) then
        allocate(q_tau(nrr_div,ndim3_tau))
        allocate(t_tau(nrr_div,ndim3_tau))
      endif
c -------- K0 on (r;tau)[q;nu] -------------------
      do ispin=1,nspin
        call k0_heg_stat(ispin,q,k0_pw(1,1,1,ispin))
      enddo
c ----------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
      do it=1,iter_vrt
        if(it>1) then
c -------- dK on (r;tau)[q;nu] -------------------
          do ispin=1,nspin
            call k_heg_stat(ispin,q,lambda_dyn(1,0,ispin),
     &                      k_pw(1,1,1,ispin),lambda_stat(1,ispin))
          enddo
        endif
        k_pw=k0_pw+k_pw
        do ispin=1,nspin
          call vertex_heg_wk_stat(k_pw(1,1,1,ispin),
     &                            lambda_dyn_tau(1,1,1,ispin),
     &                            lambda_stat(1,ispin))
        enddo
        if(spin_flips.and.it>=2) then
c -------- Get Q(q';tau)[q;v_d] ---------------------------------------
          call q_heg_stat(q_tau,k_pw)
c -------- Get T(t;v_d) -----------------------------------------
          call t_heg_stat(q,q_tau,t_tau)
          do ispin=1,nspin
            call vertex_heg_gt_stat(ispin,lambda_dyn_tau(1,1,1,ispin),
     &                              t_tau)
          enddo
        endif 
        do ispin=1,nspin
          call vertex_heg_transf_stat(lambda_dyn_tau(1,1,1,ispin),
     &                                lambda_dyn(1,0,ispin))
        enddo
        conv=.false.
        if(it==1) then
          lam_max=max(maxval(abs(lambda_dyn)),maxval(abs(lambda_stat)))
        else
          conv_1=max(maxval(abs(lambda_dyn-vrt_dyn_old)),
     &               maxval(abs(lambda_stat-vrt_stat_old)))/lam_max
          if(iq<=3) conv_table(it,iq)=max(conv_table(it,iq),conv_1)
          conv_0=conv_1
        endif
        vrt_dyn_old=lambda_dyn
        vrt_stat_old=lambda_stat
        if(it==iter_vrt) conv=.true.
        if(lam_max<1.d-13) conv=.true.
        if(conv) goto 22
      enddo   !! over it
22    deallocate(k0_pw,k_pw,lambda_dyn_tau,vrt_dyn_old,vrt_stat_old)
      if(spin_flips.and.iter_vrt>1) then
        deallocate(q_tau,t_tau)
      endif
      end