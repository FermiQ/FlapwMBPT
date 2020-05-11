      subroutine mbpt_loop
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical :: w_lda
      integer :: ind_k,k,ispin,n,ind_tau,it,i,nn2
      real*8, allocatable :: std(:,:,:),sod(:,:,:,:,:),ph(:,:),
     &                       pwng(:,:,:)
      nn2=ncmpl*nd_b_pbint(me_b+1)
! --- Decide if we use W_LDA or not --------------------------------
      w_lda=.false.
      if(w_sc_gw=='non'.and.ubi_0=='dft') w_lda=.true.
      if(w_lda.and.bse_kernel_p=='0000') w_lda=.false.
! ------------------------------------------------------------------
      if(lsclp) then
        if(.not.allocated(pw_sclp)) then
          allocate(pw_sclp(ncmpl*nrow_pw_max,ncol_pw_max,ndim3_tn,
     &                     ndim3_k(me_k+1)))
        endif
      else
        if(.not.allocated(p_w_q_tau_mm)) then
          allocate(p_w_q_tau_mm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                          ndim3_k(me_k+1)))
          allocate(p_w_q_tau_mi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                          ndim3_tn,ndim3_k(me_k+1)))
          allocate(p_w_q_tau_ii(nn2,ndim3_tn,ndim3_k(me_k+1)))
        endif
      endif
      if(.not.allocated(ws_head)) allocate(ws_head(0:n_tau/2))
      if(.not.allocated(ws_head_nu)) allocate(ws_head_nu(0:n_nu))
c -------------- G0*G0 - part --------------------------------------
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
     &                     p_w_q_tau_ii)
      endif
c          call metal_insulator(ph)
      if(.not.allocated(p_head_all)) allocate(p_head_all(6,ndim3_nu))
      if(.not.allocated(ph_dir)) allocate(ph_dir(n_phead_dir,ndim3_nu))
      call p_head_num(ph,p_head_all,ph_dir)
      deallocate(ph)
c ------------------------------------------------------------------
c      allocate(pwng(ncmpl*n_pbtot,0:3,ndim3_nu))
c      call p_wing_num_prep(pwng,p_w_q_tau_mm,p_w_q_tau_mi,p_w_q_tau_ii)
      if(.not.allocated(p_wing_all)) allocate(p_wing_all(ncmpl*3,
     &                                         n_opt_pb,ndim3_nu))
      p_wing_all=0.d0
c      call p_wing_num(pwng,p_wing_all)
c      deallocate(pwng)
c ------------------------------------------------------------------
      call screened_w
c ------------------------------------------------------------------
      sig_c_tau=0.d0
      call sigx_k_solid
      call sigc_r_tau
      if(lsclp) then
        deallocate(pw_sclp)
      else
        deallocate(p_w_q_tau_mm,p_w_q_tau_mi,p_w_q_tau_ii)
      endif
c ------------------------------------------------------------------
      v_h_bb=v_h_bb_old
      call output_static_band(v_h_bb,0)
      call output_static_band(sigx_solid_k,2)
      call mbpt_spectrum
c --------- Here we find a new spectrum ----------
      call bands_sig_0
      if(nbndf_bnd_0(2)/=0) call bands_vertex
      if(nbndf_bnd_0(2)/=0.and.mesh_k_c) then
c ---- Low energy bands for vertex corrections --------------
        call bands_vertex1
c ---- Integrals for PW-susceptibilities ----------------------------
        call psi_psi_pw_red
      endif
c ---- Output of Sigma on OMEGA ------------------------------------
      allocate(sod(nbndf,2,ndim3_omega,ndim3_k(me_k+1),nspin))
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          allocate(std(n,2,ndim3_tau))
          do ind_tau=1,ndim3_tau
            do it=1,2
              do i=1,n
                std(i,it,ind_tau)=sig_c_tau(i,i,it,ind_tau,ind_k,ispin)
              enddo
            enddo
          enddo
          call tau_to_omega_ba1(sod(1,1,1,ind_k,ispin),nbndf,std,n,n)
          deallocate(std)
        enddo
      enddo
      call output_sigma_band(sod,'sigc_w','sigc',evolt/2)
      call output_sigc_rax(sod,'SigmaC')
      deallocate(sod)
c -----------------------------------------------------------------
      call search_mu_0(chem_pot,nelec,chem_pot)
      call set_g
      if(maswrk) write(iun,*)' Chemical Potential: ',chem_pot
      call timel('****** MBPT_LOOP finished **********')
      end
