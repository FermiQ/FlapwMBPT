      subroutine gw_loop
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      character*8 :: scurr,snext,ubi_1
      logical :: w_lda,w_old,w_rec,w_hf,w_lda0,w_hf0
      integer :: it0,it,nn2
      real*8, allocatable :: ph(:,:),pwng(:,:,:)
      nn2=ncmpl*nd_b_pbint(me_b+1)
c ------------------------------------------------------------------ 
      it0=iter_dft+iter_hf
c -----------------------------------------------------------------
      if(iter>it0) then
! --- Decide if we use W_LDA or not --------------------------------
! --- Decide if we use W_HF or not --------------------------------
! --- Decide if we use W_OLD or not --------------------------------
        scurr='scf'
        if(ubi==' gw'.and.w_sc_gw=='non') scurr='non'
        if(ubi==' qp'.and.w_sc_qp=='non') scurr='non'
        if(ubi=='psi'.and.w_sc_psi=='non') scurr='non'
        if(ubi=='bsp'.and.w_sc_bsp=='non') scurr='non'
        ubi_1='scf'
        if(iter+1>it0.and.iter+1<=it0+iter_gw) ubi_1=' gw'
        it=it0+iter_gw
        if(iter+1>it.and.iter+1<=it+iter_qp) ubi_1=' qp'
        it=it+iter_qp
        if(iter+1>it.and.iter+1<=it+iter_psi) ubi_1='psi'
        it=it+iter_psi
        if(iter+1>it.and.iter+1<=it+iter_bsp) ubi_1='bsp'
        snext='scf'
        if(ubi_1==' gw'.and.w_sc_gw=='non') snext='non'
        if(ubi_1==' qp'.and.w_sc_qp=='non') snext='non'
        if(ubi_1=='psi'.and.w_sc_psi=='non') snext='non'
        if(ubi_1=='bsp'.and.w_sc_bsp=='non') snext='non'
        if(iter+1>it+iter_bsp) snext='scf'
        w_lda=.false.
        w_hf=.false.
        w_old=.false.
        if(ubi==' gw') then
          if(w_sc_gw=='non') then
            if(ubi_0=='dft') w_lda=.true.
            if(ubi_0==' hf') w_hf=.true.
            if(ubi_0==' gw') w_old=.true.
          else if(w_sc_gw=='lda') then
            w_lda=.true.
          endif
        else if(ubi==' qp') then
          if(w_sc_qp=='non') then
            if(ubi_0=='dft') w_lda=.true.
            if(ubi_0==' hf') w_hf=.true.
            if(ubi_0==' gw'.and.w_sc_gw=='non') w_old=.true.
            if(ubi_0==' qp') w_old=.true.
          else if(w_sc_qp=='lda') then
            w_lda=.true.
          endif
        else if(ubi=='psi') then
          if(w_sc_psi=='non') then
            if(ubi_0=='dft') w_lda=.true.
            if(ubi_0==' hf') w_hf=.true.
            if(ubi_0==' gw'.and.w_sc_gw=='non') w_old=.true.
            if(ubi_0==' qp'.and.w_sc_qp=='non') w_old=.true.
            if(ubi_0=='psi') w_old=.true.
          endif
        else if(ubi=='bsp') then
          if(w_sc_bsp=='non') then
            if(ubi_0=='dft') w_lda=.true.
            if(ubi_0==' hf') w_hf=.true.
            if(ubi_0==' gw'.and.w_sc_gw=='non') w_old=.true.
            if(ubi_0==' qp'.and.w_sc_qp=='non') w_old=.true.
            if(ubi_0=='psi'.and.w_sc_psi=='non') w_old=.true.
            if(ubi_0=='bsp') w_old=.true.
          endif
        endif
        w_lda0=w_lda
        if(w_lda.and.bse_kernel_p=='0000') w_lda=.false.
        w_hf0=w_hf
        if(w_hf.and.bse_kernel_p=='0000') w_hf=.false.
        p_q0_vrt=0.d0
! --------- dP and dSigma are calculated here ----------------------      
        if(w_hf0) then
          if(w_hf) then
            if(hybrid_type/=0) key_hybrid=1
            call vertex_p_hf(w_old)
            key_hybrid=0
          endif
        else if(.not.w_lda0) then
          if(ubi=='psi'.or.ubi=='bsp') then
            call vertex_p(w_old,2)
          endif
          if(ubi==' qp'.and.w_sc_qp=='non'.and.ubi_0==' gw'
     &       .and.w_sc_gw=='scf') then
            call vertex_p(w_old,2)
          endif
        endif
! ------------------------------------------------------------------
        if(lsclp) then
          if(.not.allocated(pw_sclp)) then
            allocate(pw_sclp(ncmpl*nrow_pw_max,ncol_pw_max,ndim3_tn,
     &                       ndim3_k(me_k+1)))
          endif
        else
          if(.not.allocated(p_w_q_tau_mm)) then
            allocate(p_w_q_tau_mm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                            ndim3_k(me_k+1)))
            allocate(p_w_q_tau_mi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                            ndim3_tn,ndim3_k(me_k+1)))
            allocate(p_w_q_tau_ii(nn2,ndim3_tn,ndim3_k(me_k+1)))
          endif
        endif
        if(.not.allocated(ws_head)) allocate(ws_head(0:n_tau/2))
        if(.not.allocated(ws_head_nu)) allocate(ws_head_nu(0:n_nu))
        if(w_old) then
          if(iter==restart_begin) then
            if(lsclp) then
              call write_pw_sclp(1)
            else
              call write_pw_q_tau_mm(1,p_w_q_tau_mm)
              call write_pw_q_tau_mi(1,p_w_q_tau_mi)
              call write_pw_q_tau_ii(1,p_w_q_tau_ii)
            endif
          endif
        else
          if(w_lda) then
            call response_00_dft
          else
c -------------- GG - part -----------------------------------------
            call p_r_tau_solid
          endif
c ------------------------------------------------------------------
          allocate(ph(n_k_head_0,ndim3_nu))
          if(lsclp) then
            call p_q0_num_prep_sclp(ph,pw_sclp)
          else
            call p_q0_num_prep(ph,p_w_q_tau_mm,p_w_q_tau_mi,
     &                         p_w_q_tau_ii)
          endif
c          call metal_insulator(ph)
          if(.not.allocated(p_head_all))
     &        allocate(p_head_all(6,ndim3_nu))
          if(.not.allocated(ph_dir))
     &        allocate(ph_dir(n_phead_dir,ndim3_nu))
          call p_head_num(ph,p_head_all,ph_dir)
          deallocate(ph)
c ------------------------------------------------------------------
c          allocate(pwng(ncmpl*n_pbtot,0:3,ndim3_nu))
c          call p_wing_num_prep(pwng,p_w_q_tau_mm,p_w_q_tau_mi,
c     &                         p_w_q_tau_ii)
          if(.not.allocated(p_wing_all))
     &        allocate(p_wing_all(ncmpl*3,n_opt_pb,ndim3_nu))
          p_wing_all=0.d0
c          call p_wing_num(pwng,p_wing_all)
c          deallocate(pwng)
c ------------------------------------------------------------------
          call screened_w
! ------- Decide if we have to store this W ------------------------
          w_rec=.false.
          if(iter==restart_end) then
            if(ubi==' gw'.and.w_sc_gw=='non') w_rec=.true.
            if(ubi==' qp'.and.w_sc_qp=='non') w_rec=.true.
            if(ubi=='psi'.and.w_sc_psi=='non') w_rec=.true.
            if(ubi=='bsp'.and.w_sc_bsp=='non') w_rec=.true.
          endif
          if(w_rec) then
            if(lsclp) then
              call write_pw_sclp(0)
            else
              call write_pw_q_tau_mm(0,p_w_q_tau_mm)
              call write_pw_q_tau_mi(0,p_w_q_tau_mi)
              call write_pw_q_tau_ii(0,p_w_q_tau_ii)
            endif
          endif
        endif
        sig_c_tau=0.d0
      endif
! --------- dP and dSigma are calculated here ----------------------
      if(ubi=='psi'.or.ubi=='bsp') then
        if(ubi=='psi') call vertex_sig_psi(psi_sig)
        if(ubi=='bsp') then
          if(iter_sigma_gwg<=1) call vertex_sig_psi(psi2_sig)
          if(iter_sigma_gwg>1) call vertex_sig_inf(psi2_sig)
        endif
      endif
! ------------------------------------------------------------------
c --- Here we get GW part of Sigma in TAU, mix it -----------
c ----------- and we get it on OMEGA (if need) -------      
      call self_energy
      if(iter>it0) then
        if(lsclp) then
          if(scurr=='scf') then
            deallocate(pw_sclp)
          else if(scurr=='non') then
            if(snext=='scf') then
              deallocate(pw_sclp)
            endif
          endif
          if(iter==restart_end) then
            if(allocated(pw_sclp)) deallocate(pw_sclp)
          endif
        else
          if(scurr=='scf') then
            deallocate(p_w_q_tau_mm,p_w_q_tau_mi,p_w_q_tau_ii)
          else if(scurr=='non') then
            if(snext=='scf') then
              deallocate(p_w_q_tau_mm,p_w_q_tau_mi,p_w_q_tau_ii)
            endif
          endif
          if(iter==restart_end) then
            if(allocated(p_w_q_tau_mm)) deallocate(p_w_q_tau_mm)
            if(allocated(p_w_q_tau_mi)) deallocate(p_w_q_tau_mi)
            if(allocated(p_w_q_tau_ii)) deallocate(p_w_q_tau_ii)
          endif
        endif
      endif
c ---- We get full Sigma on OMEGA ----------------------------------
      call self_energy_1
c ------------------------------------------------------------------
      if((ubi=='psi'.or.ubi=='bsp').and..not.metal) then
        call bands_gw_qp(chem_pot)
        call hovb_lucb(e_qp,hovb,lucb)
        chem_pot=0.5d0*(hovb+lucb)
      else
	    call search_mu_0(chem_pot,nelec,chem_pot)
	  endif
      if(maswrk) write(iun,*)' Chemical Potential: ',chem_pot
      if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
        call bands_gw_qp(chem_pot)
      endif
c ------------------------------------------------------------------
	  call set_g
      call timel('******** GW_LOOP finished **********')
      end
