      subroutine vertex_heg_a_spl(iter_vrt,iq,ind_nu,nom,lambda_dyn,
     &                            lambda_dyn_sig2,
     &                            lambda_stat_nu,lambda_stat_nu_sig2,
     &                            g_om,g_om_x,wexa,wspl_geom,wspl_asy,
     &                            wasy,conv_table,
     &                            conv_1,spin_flips,ntr,q,kma)
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
      integer, intent(in) :: iter_vrt,nom,ind_nu,iq,ntr,kma
      real*8, intent(in) :: wexa(0:n_nu_exa,nstar_c),q(*),
     &                      wspl_geom(0:n_nu_geom+2,4,nstar_c),
     &                      wspl_asy(n_nu_asy+1,4,nstar_c),
     &                      wasy(nstar_c)
      complex*16, intent(in) :: g_om(nstar_c,nom,2,nspin),
     &                          g_om_x(nstar_c,nom,2,nspin)
      real*8, intent(inout) :: conv_table(2:ntr,0:7,kma)
      real*8, intent(out) :: conv_1
      complex*16, intent(out) :: lambda_stat_nu(nrr_div_c,nspin),
     &                           lambda_stat_nu_sig2(nrr_div_c,nspin),
     &                           lambda_dyn(nrr_div_c,nom,nspin),
     &                           lambda_dyn_sig2(nrr_div_c,nom,nspin)
      logical :: conv
      integer :: it,i_nu,i1_nu,i,nom_nu,ispin,iqc
      real*8 :: w1,conv_0,lam_max
      real*8, allocatable :: w_nu_nu(:,:,:)
      complex*16, allocatable :: lambda_dyn_tau(:,:,:),
     &                           k0_pw(:,:,:),k_pw(:,:,:),
     &                           q_tau(:,:),t_tau(:,:),
     &                           vrt_dyn_old(:,:,:),
     &                           vrt_stat_old(:,:)
      conv_1=0.d0
      iqc=k_npnt_in_line_c(iq)
c --------------------------------------------------------------------
	i_nu=me_t*ndim3_nu+ind_nu-1
      lambda_stat_nu=(0.d0,0.d0)
      lambda_stat_nu_sig2=(0.d0,0.d0)
      lambda_dyn=(0.d0,0.d0)
      lambda_dyn_sig2=(0.d0,0.d0)
	allocate(k0_pw(nrr_div_c,0:n_tau,nspin))
	allocate(k_pw(nrr_div_c,0:n_tau,nspin))
	allocate(lambda_dyn_tau(nrr_div_c,0:n_tau,nspin))
	allocate(vrt_dyn_old(nrr_div_c,nom,nspin))
	allocate(vrt_stat_old(nrr_div_c,nspin))
c ------ W-interpolation ------------------------------------------
      if(spin_flips) then
	  allocate(w_nu_nu(nstar_c,2,n_nu2_max))
	  do i=1,nstar_c
          do i1_nu=1,num_nu2_adapt(i_nu) 
            w1=w_nu_adapt_nu(i1_nu,i_nu)
            call interp_w_nu_spl(wexa(0,i),wspl_geom(0,1,i),
     &                           wspl_asy(1,1,i),wasy(i),
     &                           w_nu_nu(i,1,i1_nu),w1,1)
            w1=w_nu_adapt_nu(i1_nu,i_nu)-w_nu(i_nu)
            call interp_w_nu_spl(wexa(0,i),wspl_geom(0,1,i),
     &                           wspl_asy(1,1,i),wasy(i),
     &                           w_nu_nu(i,2,i1_nu),w1,1)
            w_nu_nu(i,1,i1_nu)=w_nu_nu(i,1,i1_nu)+v_q_heg_c(i)
            w_nu_nu(i,2,i1_nu)=w_nu_nu(i,2,i1_nu)+v_q_heg_c(i)
          enddo
        enddo
        allocate(q_tau(nrr_div_c,0:n_tau/2))
        allocate(t_tau(nrr_div_c,0:n_tau/2))
      endif
c ---------------------------------------------------------      
	nom_nu=num_nu2_adapt(i_nu)
c -------- K0 on (r;tau)[q;nu] -------------------
      do ispin=1,nspin
        call k0_heg(ispin,q,ind_nu,nom,k0_pw(1,0,ispin),
     &              g_om(1,1,1,ispin),g_om_x(1,1,1,ispin))
      enddo
c ----------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
      do it=1,iter_vrt
        if(it>1) then
c -------- dK on (r;tau)[q;nu] -------------------
          do ispin=1,nspin
            call k_heg(ispin,q,nom,ind_nu,lambda_dyn(1,1,ispin),
     &                 k_pw(1,0,ispin),lambda_stat_nu(1,ispin),
     &                 g_om(1,1,1,ispin),g_om_x(1,1,1,ispin))
          enddo
        endif
        k_pw=k0_pw+k_pw
        do ispin=1,nspin
          call vertex_heg_wk(q,k_pw(1,0,ispin),
     &                       lambda_dyn_tau(1,0,ispin),
     &                       lambda_stat_nu(1,ispin))
        enddo
        if(it==1) then
          do ispin=1,nspin
            call vertex_heg_transf(q,ind_nu,nom,
     &                             lambda_dyn_tau(1,0,ispin),
     &                             lambda_dyn(1,1,ispin))
          enddo
          lambda_dyn_sig2=lambda_dyn
          lambda_stat_nu_sig2=lambda_stat_nu
        endif
        if(spin_flips) then
c -------- Get Q(q';tau)[q;v_d] ---------------------------------------
          call q_heg(q_tau,k_pw)
c -------- Get T(t;v_d) -----------------------------------------
          call t_heg(q,q_tau,t_tau,ind_nu,w_nu_nu)
          do ispin=1,nspin
            call vertex_heg_gt(ispin,lambda_dyn_tau(1,0,ispin),t_tau)
          enddo
        endif
        do ispin=1,nspin
          call vertex_heg_transf(q,ind_nu,nom,
     &                           lambda_dyn_tau(1,0,ispin),
     &                           lambda_dyn(1,1,ispin))
        enddo 
        conv=.false.
        if(it==1) then
          lam_max=max(maxval(abs(lambda_dyn)),
     &                maxval(abs(lambda_stat_nu)))
        else
          conv_1=max(maxval(abs(lambda_dyn-vrt_dyn_old)),
     &               maxval(abs(lambda_stat_nu-vrt_stat_old)))/lam_max
          if(i_nu<=7.and.iqc<=kma.and.iqc/=0) conv_table(it,i_nu,iqc)=
     &                            max(conv_table(it,i_nu,iqc),conv_1)
          conv_0=conv_1
        endif
        vrt_dyn_old=lambda_dyn
        vrt_stat_old=lambda_stat_nu
        if(it==iter_vrt) conv=.true.
        if(lam_max<1.d-13) conv=.true.
        if(conv) then
          lambda_dyn_sig2=lambda_dyn-lambda_dyn_sig2
          lambda_stat_nu_sig2=lambda_stat_nu-lambda_stat_nu_sig2
          goto 22
        endif
      enddo   !! over it
22    deallocate(k0_pw,k_pw,lambda_dyn_tau,vrt_dyn_old,vrt_stat_old)
      if(spin_flips) then
        deallocate(w_nu_nu,q_tau,t_tau)
      endif
      end