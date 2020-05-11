      subroutine memory1_solid
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer :: me1,nlmb
      me1=me+1
      nlmb=nrel*maxb+1
      allocate(g_loc_0(maxel,maxel,natom,nspin_0))
      allocate(ro(maxmt))
      allocate(p_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      allocate(q_f(0:maxnrad,maxntle,nlmb,nsort,nspin))
      allocate(gfun_old(maxwf,nspin))
      allocate(gfund_old(maxwf,nspin))
      allocate(t_ro(maxlfun,maxlfun,maxnsym,nsort,nspin_0))
      allocate(v_r0_asym(0:2,maxnsym,nspin_0,nsort))
      allocate(ro_r0_asym(0:2,maxnsym,nspin_0,nsort))
      allocate(z_r0_asym(0:2,nsort))
      allocate(v_mt_h(maxmt/nspin_0))
      allocate(v_mt_h_old(maxmt/nspin_0))
      allocate(v_mt_xc(maxmt))
      v_mt_xc=0.d0
      allocate(pv(maxel,maxel,natom,nspin_0))
      allocate(pvj(maxel,maxel,natom))
      allocate(g_full_00(nbndf,ndim3_k(me_k+1),nspin_0))
      if(iter_gw+iter_qp+iter_psi+iter_bsp/=0) then
        allocate(g_full(nbndf,nbndf,2,ndim3_tau,ndim3_k(me_k+1),
     &                  nspin_0))
        allocate(sig_c_tau(nbndf,nbndf,2,ndim3_tau,ndim3_k(me_k+1),
     &                     nspin_0))
        allocate(sig_c_tau_old(nbndf,nbndf,2,ndim3_tau,
     &                         ndim3_k(me_k+1),nspin_0))
        sig_c_tau_old=0.d0
      endif
      if(iter_max_scf>iter_dft+iter_hf) then
        allocate(e_qp(nbndf,npnt,nspin_0))
        if(iter_qp/=0) then
          allocate(z_qp(nfun,nbndf,ndim3_k(me_k+1),nspin_0))
          allocate(q_qp(nbndf,nbndf,ndim3_k(me_k+1),nspin_0))
          allocate(qp_spectr_weight(nbndf,ndim3_k(me_k+1),nspin_0))
        else
          allocate(z_qp(nfun,min(nbndf,9),ndim3_k(me_k+1),nspin_0))
        endif
      endif
      if(iter_max_scf>iter_dft) then
        allocate(v_h_bb_old(nd_tb_ij(me_t_b+1),ndim3_k(me_k+1),nspin_0))
        allocate(v_xc_bb(nd_tb_ij(me_t_b+1),ndim3_k(me_k+1),nspin_0))
        allocate(v_xc_bb_old(nd_tb_ij(me_t_b+1),ndim3_k(me_k+1),
     &                       nspin_0))
        allocate(v_h_bb(nd_tb_ij(me_t_b+1),ndim3_k(me_k+1),nspin_0))
        allocate(sigx_solid_k(nd_tb_ij(me_t_b+1),ndim3_k(me_k+1),
     &                        nspin_0))
        allocate(sigx_solid_k_old(nd_tb_ij(me_t_b+1),ndim3_k(me_k+1),
     &                            nspin_0))
        sigx_solid_k_old=0.d0
      endif
      if(magn.eq.2) then
        allocate(spmt(maxmtb))
        allocate(t_m(maxlfun,maxlfun,2,maxnsymb,nsort))
        spmt=0.d0
      endif
      allocate(ntle_pb(0:maxpb,nsort))
      end
