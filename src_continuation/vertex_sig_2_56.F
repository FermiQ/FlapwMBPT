      subroutine vertex_sig_2_56(sig_56_2,psi2)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*3, intent(in) :: psi2
      complex*16, intent(out) :: sig_56_2(nbndf_bnd,nbndf_bnd,2,
     &                                    ndim3_tau,ndim_k_red(me_k+1),
     &                                    nspin)
      integer :: k,ispin,iq,ind_nu,i_nu,i,nom,s2,i_omega,np,mode,
     &           ind_q,nbb,ib,ind_omega,mx2,i_tau,kq,kst,iq_line,it,
     &           nomnu,ind_s2,kf,iqf,ka,j
      real*8 :: de,green0,v(3),gtild(3,2)
      complex*16 :: cdum
      real*8, allocatable :: gexa(:,:,:,:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:,:,:),
     &                       gasy(:,:,:,:,:),g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),fif(:,:,:,:,:),
     &                       fif_full(:,:,:,:,:)
      complex*16, allocatable :: k0_pw(:,:,:,:,:,:),trn(:,:),
     &                           g_om_nu(:,:,:,:,:,:),v2_q(:,:,:),
     &                           z_red(:,:,:,:),ps3(:,:,:,:,:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           s_56(:,:,:,:,:),aw56(:,:,:,:,:),
     &                           phf(:,:),
     &                           z_full(:,:,:,:),a_full(:,:,:,:)
      mx2=maxel_red**2
      mode=3
      if(psi2=='dyn') mode=4
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
      allocate(fif_full(nrel,maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ---------- Interaction matrices in Q-orthonormal product basis ----      
      allocate(v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c))
      v2_q=(0.d0,0.d0)
      if(psi2=='sta') then
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          if(i_nu/=nu_w_stat) cycle
          do k=1,npnt_c
            ka=k_a_from_c(k)
            do kst=1,k_star_c(k)
              iq=k_list_c(kst,k)
              v2_q(:,:,iq)=v_red_q(:,:,ka)+w_red_q(:,:,ind_nu,ka)
              call sym_w_red(iq,k,v2_q(1,1,iq))
            enddo   !! over kst
          enddo   !! over k
        enddo
        if(nproc_t/=1) call dgop(v2_q,2*n_pbtot_red**2*nqdiv_c,'  +',
     &                             comm_t)
      else if(psi2=='dyn') then
        call set_wnu_vrt
c ---------------------------------------------------------------------
        do k=1,npnt_c
          ka=k_a_from_c(k)
          do kst=1,k_star_c(k)
            iq=k_list_c(kst,k)
            v2_q(:,:,iq)=v_red_q(:,:,ka)
            call sym_w_red(iq,k,v2_q(1,1,iq))
          enddo   !! over kst
        enddo   !! over k
      endif
      call timel('*** VERTEX_BND_GW : all W finished *')
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
c -----------------
      allocate(g_rs_mm(nrel,nfun_red,nfun_red,0:n_tau,nqdiv_c,nspin))
      allocate(g_rs_mi(nrel,nfun_red,nrel*nr_full_red,0:n_tau,nqdiv_c,
     &                 nspin))
      allocate(g_rs_im(nrel,nrel*nr_full_red,nfun_red,0:n_tau,nqdiv_c,
     &                 nspin))
      allocate(g_rs_ii(nrel,nrel*nr_red_red,nrel*nr_full_red,0:n_tau,
     &                 nqdiv_c,nspin))
      do ispin=1,nspin
        call g_real_space(g_tau(1,1,1,1,1,ispin),
     &                    g_rs_mm(1,1,1,0,1,ispin),
     &                    g_rs_mi(1,1,1,0,1,ispin),
     &                    g_rs_im(1,1,1,0,1,ispin),
     &                    g_rs_ii(1,1,1,0,1,ispin),z_red(1,1,1,ispin),
     &                    a_red(1,1,1,ispin),phf,ispin)
      enddo
      call timel('**** G_REAL_SPACE finished *********')
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
      allocate(k0_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2,nqdiv_c,nspin))
      allocate(vdy_store(nrs_red,nrs_red,0:n_tau,2,nqdiv_c,nspin))
      allocate(eqr(nqdiv_c))
c --------------------------------------------------------------------
      allocate(aw56(nra_red,nra_red,nqdiv_c,0:n_tau,nspin))
      aw56=(0.d0,0.d0)
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
	    call eqr_vertex_line(eqr,iq)
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
            call ppm_factors_full(k,kq,iq,fif_full(1,1,1,1,ispin),
     &	                          z_full(1,1,1,ispin),
     &                            a_full(1,1,1,ispin),
     &                            ps3(1,1,1,k,ispin),gtild,
     &                            n_low_bnd(1,ispin))
          enddo   !! over k
        enddo   !! over ispin
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          nom=num_omega1_adapt(i_nu)
          nomnu=num_nu2_adapt(i_nu)
          allocate(g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,npnt_c,nspin))
          g_om_nu=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
          call g_interp_bnd(gexa,gspl_geom,gspl_asy,gasy,g_om_nu,i_nu,
     &                      nom)
c --------------------------------------------------------------------
          do ind_s2=1,ndim4_pbr(me_b+1,iq)
            s2=n4_mpi_pbr(me_b+1,iq)+ind_s2
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND starts *****************')
c -------- Get zero order generalized susceptibility -------------------
            call k0_bnd(ind_nu,iq,k0_pw,g_om_nu,ps3,s2,g_x,ex,nom)
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND finished ***************')
c ----------------------------------------------------------------------
            vdy_store=(0.d0,0.d0)
            call vertex_bnd_loop(.true.,iq,ind_nu,1,2,s2,nom,
     &                           nomnu,k0_pw,fif0,z_red,a_red,phf,ax,
     &                           cdum,cdum,cdum,cdum,cdum,cdum,
     &                           cdum,cdum,v2_q,cdum,cdum,nspin,mode)
c --------------------------------------------------------------------
            do ispin=1,nspin
              call self_energy_1_dyn_a(ispin,iq,ind_nu,nom,
     &                                 g_om_nu(1,1,1,1,1,ispin),
     &                                 v_red_q(1,1,ka),s2,
     &                                 w_red_q(1,1,ind_nu,ka),
     &                                 vdy_store(1,1,0,1,1,ispin),
     &                                 aw56(1,1,1,0,ispin),
     &                                 z_red(1,1,1,ispin),
     &                                 a_red(1,1,1,ispin),
     &                                 fif(1,1,1,1,ispin),phf,psi2)
            enddo
            if(iq==1.and.i_nu==0.and.s2==1)
     &          call timel('**** Self En_1_dyn_A finished ******')
          enddo   !! over ind_s2
          if(iq==1.and.i_nu==0) 
     &       call timel('**** S2-loop finished **************')
          deallocate(g_om_nu)
        enddo   !! over ind_nu
        if(iq==1) call timel('**** Nu-loop finished **************')
        if(iq==1) call timel('**** END of Q=1 loop ***************')
      enddo  !! over ind_q
      call timel('**** END of Q loop *****************')
      deallocate(v2_q,k0_pw,vdy_store,gexa,gspl_geom,gspl_asy,eqr,
     &           gasy,fif,g_x,ex,fif_full,z_full,a_full,g_rs_mm,g_rs_mi,
     &           g_rs_im,g_rs_ii)
      if(psi2=='dyn') deallocate(wnu_vrt)
c ---------------------------------------------------------------      
      deallocate(g_tau)
      if(nproc_k/=1) call dgop(aw56,
     &                         2*nra_red**2*(n_tau+1)*nqdiv_c*nspin,
     &                         '  +',comm_k)
      if(nproc_b/=1) call dgop(aw56,
     &                           2*nra_red**2*(n_tau+1)*nqdiv_c*nspin,
     &                           '  +',comm_b)
      if(nproc_t/=1) call dgop(aw56,
     &                           2*nra_red**2*(n_tau+1)*nqdiv_c*nspin,
     &                           '  +',comm_t)
      allocate(s_56(nra_red,nra_red,nqdiv_c,0:n_tau,nspin))
      s_56=(0.d0,0.d0)
      do ispin=1,nspin
        call self_energy_1_dyn_b(aw56(1,1,1,0,ispin),
     &                           s_56(1,1,1,0,ispin))
      enddo
      deallocate(aw56)
c ------- Symmetrization of Sigma + -----> band repr.  ------------
      call sym_sig_k_red_a(sig_56_2,s_56,z_red,a_red)
      deallocate(s_56,phf,nind_fif,ind_fif,fif0,ax)
      call timel('**** VERTEX_SIG_2_56 finished ******')
      end
