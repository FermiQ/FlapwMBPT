      subroutine vertex_sig_2_14(sig_1_2)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(out) :: sig_1_2(nbndf_bnd,nbndf_bnd,2,
     &                                   ndim3_tau,ndim_k_red(me_k+1),
     &                                   nspin)
      integer :: k,ispin,iq,ind_nu,i_nu,i,nom,it_vrt,s2,i_omega,np,
     &           ind_q,nbb,ib,ind_omega,mx2,i_tau,kq,kst,iq_line,it,
     &           ind_tau,nomnu,ind_s2,kf,iqf,ka,j,nrr3
      real*8 :: de,green0,v(3),gtild(3,2)
      integer, allocatable :: indrr3(:,:,:)
      real*8, allocatable :: v4_r(:,:,:),w4_r(:,:,:,:),
     &                       gexa(:,:,:,:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:,:,:),
     &                       gasy(:,:,:,:,:),fif(:,:,:,:,:),
     &                       w_tau(:,:,:,:),g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),vmi(:,:,:),
     &                       wmi(:,:,:,:),vim(:,:,:),
     &                       wim(:,:,:,:),vii(:),
     &                       wii(:,:),fif_full(:,:,:,:)
      complex*16, allocatable :: trn(:,:),k_pw(:,:,:,:,:,:),
     &                           lambda_dyn(:,:,:,:,:,:),
     &                           lambda_stat(:,:,:,:),
     &                           g_om_nu(:,:,:,:,:,:),v2_q(:,:,:),
     &                           z_red(:,:,:,:),ps3(:,:,:,:,:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           k_rs(:,:,:,:,:,:),
     &                           s_12(:,:,:,:,:),aw12(:,:,:,:,:),
     &                           phf(:,:),
     &                           vrt2_tau(:,:,:,:,:),
     &                           z_full(:,:,:,:),a_full(:,:,:,:),
     &                           vrt2_nu(:,:,:,:)
      mx2=maxel_red**2
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
c ------------------------------------------------------------------
      allocate(fif(nrel,mx2,n_pbmtm_red,nsort,nspin))
      call aux_vrt_add(fif)
c ------------------------------------------------------------------
      call norma_wf_red
c ------ Here for zero-order vertex in full basis ------------------      
      allocate(fif_full(maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ---------- Interaction matrices in Q-orthonormal product basis ----
      allocate(w_tau(n_pbtot_red,n_pbtot_red,ndim3_tau,nqdiv_c))
      call w_tau_prepare_red(w_tau)
c ---------------------------------------------------------------------      
      allocate(v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c))
      v2_q=(0.d0,0.d0)
      do k=1,npnt_c
        ka=k_a_from_c(k)
        do kst=1,k_star_c(k)
          iq=k_list_c(kst,k)
          v2_q(:,:,iq)=v_red_q(:,:,ka)
          call sym_w_red(iq,k,v2_q(1,1,iq))
        enddo   !! over kst
      enddo   !! over k
      call timel('**** VERTEX_BND_GW : v2_q finished *')
c --------- Interaction matrices in 2-site (MM) form ------------------
      allocate(v4_r(n_pbmt_red,n_pbmt_red,nqdiv_c))
      allocate(w4_r(n_pbmt_red,n_pbmt_red,nqdiv_c,0:n_tau/2))
      call w4_for_bnd_mm(v4_r,w4_r,w_tau,v2_q,n_pbtot_red)
c --------- Interaction matrices in (MI) form ------------------
      allocate(vmi(n_pbmt_red,nr_full_red,nqdiv_c))
      allocate(wmi(n_pbmt_red,nr_full_red,nqdiv_c,0:n_tau/2))
      call w4_for_bnd_mi(vmi,wmi,w_tau,v2_q,phf)
c --------- Interaction matrices in (IM) form ------------------
      allocate(vim(nr_full_red,n_pbmt_red,nqdiv_c))
      allocate(wim(nr_full_red,n_pbmt_red,nqdiv_c,0:n_tau/2))
      call w4_for_bnd_im(vim,wim,w_tau,v2_q,phf)
c --------- Interaction matrices in (II) form ------------------
      nrr3=nr_red_red*nr_full_red*nqdiv_c
      allocate(indrr3(nr_full_red,nr_full_red,nqdiv_c))
      allocate(vii(nrr3))
      allocate(wii(nrr3,0:n_tau/2))
      call w4_for_bnd_ii(vii,wii,w_tau,v2_q,phf,nrr3,indrr3)
      deallocate(w_tau)
      call timel('*** VERTEX_BND_GW : all W finished *')
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
c ---- Exchange eigen values relative to Chem_POT ------------------
      allocate(ex(nbndf_bnd,npnt_c,nspin))
      ex=0.d0
      do ispin=1,nspin
        do k=1,npnt_c
          kf=k_a_from_c(k)
          do i=1,n_low_bnd(kf,ispin)
            ib=ind_bands_bnd(i,kf,ispin)
            ex(i,k,ispin)=e_bnd(ib,kf,ispin)-chem_pot
          enddo
        enddo   !! over k
      enddo  !! over ispin      
c --------- Exchange Green's function ------------------------------
      allocate(g_x(nbndf_bnd,0:n_tau,npnt_c,nspin))
      g_x=0.d0
      do ispin=1,nspin
        do k=1,npnt_c
          kf=k_a_from_c(k)
          do i=1,n_low_bnd(kf,ispin)
            de=ex(i,k,ispin)
            do i_tau=0,n_tau
              g_x(i,i_tau,k,ispin)=green0(de,tau_mesh(i_tau))
            enddo
          enddo
        enddo   !! over k
      enddo  !! over ispin
c -------------	
      allocate(gexa(nbndf_bnd,nbndf_bnd,2,0:n_omega_exa,npnt_c,nspin))
      allocate(gspl_geom(nbndf_bnd,nbndf_bnd,2,0:n_omega_geom+2,4,
     &                   npnt_c,nspin))
      allocate(gspl_asy(nbndf_bnd,nbndf_bnd,2,n_omega_asy+1,4,npnt_c,
     &                  nspin))
      allocate(gasy(nbndf_bnd,nbndf_bnd,2,npnt_c,nspin))
      gexa=0.d0
      gspl_geom=0.d0
      gspl_asy=0.d0
      gasy=0.d0
c --------- For Interpolation of G(w') ---------------------------------
      nbb=nbndf_bnd**2
      allocate(g0(0:n_omega,2))
      do ispin=1,nspin
        do k=1,npnt_c
          do j=1,nbndf_bnd
            do i=1,nbndf_bnd
              g0=0.d0
              do ind_omega=1,ndim3_omega
                i_omega=me_t*ndim3_omega+ind_omega-1
                do it=1,2
                  g0(i_omega,it)=g_omega(i,j,it,ind_omega,k,ispin)
                enddo
              enddo
              if(nproc_t/=1) call dgop(g0,2*(n_omega+1),'  +',
     &                                   comm_t)
              do it=1,2
                call spline_inhmg(w_omega(n_omega_exa-1),
     &                            g0(n_omega_exa-1:,it),
     &                            gspl_geom(i,j,it,:,1,k,ispin),
     &                            gspl_geom(i,j,it,:,2,k,ispin),
     &                            gspl_geom(i,j,it,:,3,k,ispin),
     &                            gspl_geom(i,j,it,:,4,k,ispin),
     &                            n_omega_geom+3,0,0.d0,0.d0)
                call spline_inhmg(xm_omega(0),
     &                            g0(n_omega_exa+n_omega_geom:,it),
     &                            gspl_asy(i,j,it,:,1,k,ispin),
     &                            gspl_asy(i,j,it,:,2,k,ispin),
     &                            gspl_asy(i,j,it,:,3,k,ispin),
     &                            gspl_asy(i,j,it,:,4,k,ispin),
     &                            n_omega_asy+1,0,0.d0,0.d0)
                gasy(i,j,it,k,ispin)=g0(n_omega,it)
                gexa(i,j,it,:,k,ispin)=g0(0:n_omega_exa,it)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(g_omega,g0)
      call timel('*** VERTEX_BND_GW : G_0 finished ***')
      allocate(k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2,nqdiv_c,nspin))
      allocate(k_rs(nrs_red,nrs_red,0:n_tau,2,nqdiv_c,nspin))
      allocate(vdy_store(nrs_red,nrs_red,0:n_tau,2,nqdiv_c,nspin))
      allocate(lambda_stat(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
c --------------------------------------------------------------------
      allocate(aw12(nra_red,nra_red,nqdiv_c,0:n_tau,nspin))
      aw12=(0.d0,0.d0)
      allocate(ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,nqdiv_c,nspin))
      call timel('*** most allocations finished ******')
c ------------------------------------------------------------------ 
      allocate(trn(n_omega1_max,0:n_nu))
      trn=(0.d0,0.d0)
      if(me_t==0) then
        do i_nu=0,n_nu
          trn(:,i_nu)=tau_from_omega_nu(:,1,1,i_nu)
        enddo
      endif
      if(nproc_t/=1) call dgop(trn,2*n_omega1_max*(n_nu+1),'  +',
     &                           comm_t)
c ------------------------------------------------------------------
      call timel('**** VERTEX_BND_GW : loop starts ***')
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        iqf=k_a_from_c(iq)
        ka=k_a_from_c(iq)
        iq_line=k_npnt_in_line_c(iq)
        np=n_pbmt_red+nplwgw_red(ka)
        ps3=(0.d0,0.d0)
        do ispin=1,nspin
          do k=1,nqdiv_c
            v=pnt_c(:,k)-pnt_c(:,iq)
            call zone1_number(v,rb0_c,ndiv_c,kq)
            kq=index_k1_c(kq)    !! for K-Q
            gtild(:,1)=0.d0
            gtild(:,2)=v-pnt_c(:,kq)
            call ppm_factors_full(k,kq,iq,fif_full(1,1,1,ispin),
     &	                          z_full(1,1,1,ispin),
     &                            a_full(1,1,1,ispin),
     &                            ps3(1,1,1,k,ispin),gtild,
     &                            n_low_bnd(1,ispin))
          enddo   !! over k
        enddo   !! over ispin
        allocate(vst_store(nrs_red,nrs_red,nqdiv_c,ndim3_nu,
     &                     ndim4_pbr(me_b+1,iq),nspin))
        vst_store=(0.d0,0.d0)
        do ind_nu=1,ndim3_nu
          ind_nu_st=ind_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          nom=num_omega1_adapt(i_nu)
          nomnu=num_nu2_adapt(i_nu)
          allocate(lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,nqdiv_c,nspin))
          allocate(g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,npnt_c,nspin))
          g_om_nu=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
          call g_interp_bnd(gexa,gspl_geom,gspl_asy,gasy,g_om_nu,i_nu,
     &                      nom)
c --------------------------------------------------------------------
          do ind_s2=1,ndim4_pbr(me_b+1,iq)
            ind_pbr_st=ind_s2
            s2=n4_mpi_pbr(me_b+1,iq)+ind_s2
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND starts *****************')
            lambda_dyn=(0.d0,0.d0)
c -------- Get zero order generalized susceptibility -------------------
            call k0_bnd(ind_nu,iq,k_pw,g_om_nu,ps3,s2,g_x,ex,nom)
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND finished ***************')
c ----------------------------------------------------------------------
            do it_vrt=1,2
              vdy_store=(0.d0,0.d0)
              call vertex_bnd_loop(.false.,iq,ind_nu,it_vrt,2,s2,nom,
     &                             nomnu,k_pw,fif0,z_red,a_red,phf,ax,
     &                             v4_r,vmi,vim,vii,w4_r,wmi,wim,wii,
     &                             v2_q,lambda_stat,lambda_dyn,nspin,2,
     &                             nrr3,indrr3)
              if(it_vrt==1) then
                do ispin=1,nspin
                  k_pw(:,:,:,:,:,ispin)=(0.d0,0.d0)
                  do k=1,nqdiv_c
c -------- Get correction to the generalized susceptibility --------
                    call k_bnd(ispin,nom,ind_nu,k,iq,
     &                         lambda_dyn(1,1,1,1,k,ispin),
     &                         k_pw(1,1,0,1,k,ispin),
     &                         lambda_stat(1,1,k,ispin),
     &                         g_om_nu(1,1,1,1,1,ispin),
     &                         g_x(1,0,1,ispin),ex(1,1,ispin))
                  enddo  !! over k
                enddo  !! over ispin
              endif
            enddo  !! over it_vrt
c -------- For dynamical self energy ----------------------------
            do ispin=1,nspin
              call self_energy_1_dyn_a(ispin,iq,ind_nu,nom,
     &                                 g_om_nu(1,1,1,1,1,ispin),
     &                                 v_red_q(1,1,iqf),s2,
     &                                 w_red_q(1,1,ind_nu,iqf),
     &                                 vdy_store(1,1,0,1,1,ispin),
     &                                 aw12(1,1,1,0,ispin),
     &                                 z_red(1,1,1,ispin),
     &                                 fif(1,1,1,1,ispin),phf,'dyn')
            enddo
            if(iq==1.and.i_nu==0.and.s2==1)
     &          call timel('**** Self En_1_dyn_A finished ******')
          enddo   !! over ind_s2
          if(iq==1.and.i_nu==0) 
     &       call timel('**** S2-loop finished **************')
          deallocate(lambda_dyn,g_om_nu)
        enddo   !! over ind_nu
        if(iq==1) call timel('**** Nu-loop finished **************')
        allocate(vrt2_tau(nbndf_bnd,nrs_red,nqdiv_c,2,ndim3_tau))
        allocate(vrt2_nu(nbndf_bnd,nrs_red,nqdiv_c,ndim3_nu))
        do ispin=1,nspin
          do ind_s2=1,ndim4_pbr(me_b+1,iq)
            s2=n4_mpi_pbr(me_b+1,iq)+ind_s2
            vrt2_tau=(0.d0,0.d0)
            call vx2_stat_transf(vst_store(1,1,1,1,ind_s2,ispin),
     &                           vrt2_tau,vrt2_nu,
     &                           z_red(1,1,1,ispin),
     &                           a_red(1,1,1,ispin),phf,ispin)
            if(iq==1.and.s2==1)
     &          call timel('**** VX2_STAT_TRANSF finished ******')
            do ind_tau=1,ndim3_tau
c --------- Static self energy of the first order in TAU --------------
              call self_energy_1_stat(ispin,ind_tau,iq,
     &                                g_tau(1,1,1,1,1,ispin),s2,
     &                                v2_q(1,1,iq),
     &                                vrt2_tau(1,1,1,1,ind_tau),
     &                                aw12(1,1,1,0,ispin),
     &                                fif(1,1,1,1,ispin),
     &                                z_red(1,1,1,ispin),
     &                                a_red(1,1,1,ispin),phf)
            enddo   !! over ind_tau
            if(iq==1.and.s2==1)
     &          call timel('**** SELF_ENERGY_1_STAT finished ***')
c --------- Semidynamic self energy of the first order ----------------
            call self_energy_1_semi(ispin,iq,g_tau(1,1,1,1,1,ispin),
     &                              s2,vrt2_nu,aw12(1,1,1,0,ispin),
     &                              fif(1,1,1,1,ispin),
     &                              z_red(1,1,1,ispin),
     &                              a_red(1,1,1,ispin),phf,
     &                              w_red_q(1,1,1,iqf))
            if(iq==1.and.s2==1)
     &          call timel('**** SELF_ENERGY_1_SEMI finished ***')
          enddo  !! over ind_s2
        enddo   !! over ispin
        deallocate(vrt2_nu,vrt2_tau,vst_store)
        if(iq==1) call timel('**** END of Q=1 loop ***************')
      enddo  !! over ind_q
      call timel('**** END of Q loop *****************')
      deallocate(v4_r,w4_r,v2_q,vmi,wmi,k_pw,vdy_store,
     &           lambda_stat,gexa,gspl_geom,gspl_asy,gasy,fif,g_x,
     &           ex,vim,wim,vii,wii,k_rs,fif_full,z_full,a_full,indrr3)
c ---------------------------------------------------------------      
      deallocate(g_tau)
      if(nproc_k/=1) call dgop(aw12,
     &                         2*nra_red**2*(n_tau+1)*nqdiv_c*nspin,
     &                         '  +',comm_k)
      if(nproc_b/=1) call dgop(aw12,
     &                           2*nra_red**2*(n_tau+1)*nqdiv_c*nspin,
     &                           '  +',comm_b)
      if(nproc_t/=1) call dgop(aw12,
     &                           2*nra_red**2*(n_tau+1)*nqdiv_c*nspin,
     &                           '  +',comm_t)
      allocate(s_12(nra_red,nra_red,nqdiv_c,0:n_tau,nspin))
      s_12=(0.d0,0.d0)
      do ispin=1,nspin
        call self_energy_1_dyn_b(aw12(1,1,1,0,ispin),
     &                           s_12(1,1,1,0,ispin))
      enddo
      deallocate(aw12)
c ------- Symmetrization of Sigma + -----> band repr.  ------------
      call sym_sig_k_red_a(sig_1_2,s_12,z_red,a_red)
      deallocate(s_12,phf,nind_fif,ind_fif,fif0,ax)
      call timel('**** VERTEX_SIG_2__14 finished *****')
      end
