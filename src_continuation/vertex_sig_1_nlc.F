      subroutine vertex_sig_1_nlc(sig_1_1,iter_s)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iter_s
      complex*16, intent(out) :: sig_1_1(nbndf_bnd,nbndf_bnd,2,
     &                                   ndim3_tau,ndim_k_red(me_k+1),
     &                                   nspin)
      integer :: k,ispin,iq,i_nu,i,nom,s2,i_omega,nbb,ib,
     &           ind_omega,mx2,i_tau,kst,it,ind_tau,n,nomnu,ind_s2,
     &           ind_k,kf,ka,j,nb,k1,k1f,k1k,i1_omega,ind_nu,
     &           ir,numtn,itr
      real*8 :: check_hermit,w1,de,green0,v(3),gtild(3),pib,
     &          herm,hm,pi2,phase
      real*8, allocatable :: grsmm(:,:,:,:,:),wnu(:,:,:,:),wexa(:,:,:),
     &                       grsmi(:,:,:,:,:),grsim(:,:,:,:,:),
     &                       grsii(:,:,:,:,:),
     &                       wspl_geom(:,:,:,:),wspl_asy(:,:,:,:),
     &                       wasy(:,:),gexa(:,:,:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:,:),
     &                       gasy(:,:,:,:),sigt(:,:,:,:),
     &                       g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),
     &                       fif_full(:,:,:,:,:),
     &                       sig_dyn_om(:,:,:,:,:,:)
      complex*16, allocatable :: k_pw(:,:,:,:,:),ax(:,:,:,:),
     &                           t_pw(:,:,:,:,:),wnuq(:,:,:,:),
     &                           g_om_om(:,:,:,:,:),
     &                           g_nu_om(:,:,:,:,:),v2_q(:,:,:),
     &                           z_red(:,:,:,:),ps3(:,:,:,:),
     &                           a_red(:,:,:,:),
     &                           t_rr_k(:,:),
     &                           z_full(:,:,:,:),a_full(:,:,:,:),
     &                           sig_om(:,:,:),phf(:,:),
     &                           t_stat(:,:,:,:,:)
      pi2=pi+pi
      pib=pi/betta_t
      mx2=maxel_red**2
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
c ---------------------------------------------------------------------      
      call norma_wf_red
c ------ Here for zero-order vertex in full basis ------------------      
      allocate(fif_full(nrel,maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ---------- Interaction matrices in Q-orthonormal product basis ----
      allocate(wnu(n_pbtot_red,n_pbtot_red,0:n_nu,npnt_c))
      wnu=0.d0
      do iq=1,npnt_c
        ka=k_a_from_c(iq)
        nb=n_pbmt_red+nplwgw_red(ka)
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
	      call pack_hermit(w_red_q(1,1,ind_nu,ka),wnu(1,1,i_nu,iq),nb,
     &                     n_pbtot_red,n_pbtot_red,0.d0,1.d0)
        enddo  !! over ind_nu
      enddo  !! over iq
      if(nproc_t/=1) call dgop(wnu,n_pbtot_red**2*(n_nu+1)*npnt_c,
     &                           '  +',comm_t)
c ------ W-interpolation - on regullar NU for Semi-Sigma ---------
      allocate(wexa(n_pbtot_red,n_pbtot_red,0:n_nu_exa))
      allocate(wspl_geom(n_pbtot_red,n_pbtot_red,0:n_nu_geom+2,4))
      allocate(wspl_asy(n_pbtot_red,n_pbtot_red,n_nu_asy+1,4))
      allocate(wasy(n_pbtot_red,n_pbtot_red))
      allocate(wnuq(n_pbtot_red,n_pbtot_red,ndim3_nu,nqdiv_c))
      do iq=1,nqdiv_c
        call w_spl_q_red(iq,wnu,wexa,wspl_geom,wspl_asy,wasy)
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          w1=w_nu(i_nu)
          call interp_w_nu_sq_spl(wexa,wspl_geom,wspl_asy,wasy,
     &                            wnuq(1,1,ind_nu,iq),w1,n_pbtot_red)
        enddo
      enddo
      deallocate(wexa,wspl_geom,wspl_asy,wasy)
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
      call timel('**** VERTEX_SIG_1 : v2_q finished **')
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
c -----------------
      allocate(grsmm(nrel,nfun_red,nfun_red,0:n_tau,nqdiv_c))
      allocate(grsmi(nrel,nfun_red,nrel*nr_full_red,0:n_tau,nqdiv_c))
      allocate(grsim(nrel,nrel*nr_full_red,nfun_red,0:n_tau,nqdiv_c))
      allocate(grsii(nrel,nrel*nr_red_red,nrel*nr_full_red,0:n_tau,
     &               nqdiv_c))
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
      allocate(t_rr_k(nqdiv_c,nqdiv_c))
      do k=1,nqdiv_c
	    do ir=1,nqdiv_c
	      phase=-pi2*dot_product(pnt_c(:,k),rr_coarse(:,ir))
	      t_rr_k(ir,k)=dcmplx(cos(phase),sin(phase))
	    enddo
	  enddo
c -------------	
      allocate(gexa(nbndf_bnd,nbndf_bnd,2,0:n_omega_exa,npnt_c))
      allocate(gspl_geom(nbndf_bnd,nbndf_bnd,2,0:n_omega_geom+2,4,
     &                   npnt_c))
      allocate(gspl_asy(nbndf_bnd,nbndf_bnd,2,n_omega_asy+1,4,npnt_c))
      allocate(gasy(nbndf_bnd,nbndf_bnd,2,npnt_c))
      call timel('*** VERTEX_SIG_1 : G_0 finished ****')
      allocate(k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c,0:n_tau,2))
c --------------------------------------------------------------------
      allocate(ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,nqdiv_c))
      allocate(sig_dyn_om(nbndf_bnd,nbndf_bnd,2,ndim3_omega,
     &                    ndim_k_red(me_k+1),nspin))
      sig_dyn_om=0.d0
      sig_1_1=(0.d0,0.d0)
      call timel('*** most allocations finished ******')
c ------------------------------------------------------------------
      do ispin=1,nspin
        call g_real_space(g_tau(1,1,1,1,1,ispin),grsmm,grsmi,
     &                    grsim,grsii,z_red(1,1,1,ispin),
     &                    a_red(1,1,1,ispin),phf,ispin)
        call timel('**** G_REAL_SPACE finished *********')
        gexa=0.d0
        gspl_geom=0.d0
        gspl_asy=0.d0
        gasy=0.d0
c --------- For Interpolation of G(w') ---------------------------------
        nbb=nbndf_bnd**2
        allocate(g0(0:n_omega,2))
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
     &                            gspl_geom(i,j,it,:,1,k),
     &                            gspl_geom(i,j,it,:,2,k),
     &                            gspl_geom(i,j,it,:,3,k),
     &                            gspl_geom(i,j,it,:,4,k),
     &                            n_omega_geom+3,0,0.d0,0.d0)
                call spline_inhmg(xm_omega(0),
     &                            g0(n_omega_exa+n_omega_geom:,it),
     &                            gspl_asy(i,j,it,:,1,k),
     &                            gspl_asy(i,j,it,:,2,k),
     &                            gspl_asy(i,j,it,:,3,k),
     &                            gspl_asy(i,j,it,:,4,k),
     &                            n_omega_asy+1,0,0.d0,0.d0)
                gasy(i,j,it,k)=g0(n_omega,it)
                gexa(i,j,it,:,k)=g0(0:n_omega_exa,it)
              enddo
            enddo
          enddo
        enddo
        deallocate(g0)
        call timel('**** VERTEX_SIG_1 : loop starts ****')
        do ind_k=1,ndim_k_red(me_k+1)
          k1=n_mpi_k_red(me_k+1)+ind_k
          k1f=k_a_from_c(k1)
          nb=n_low_bnd(k1f,ispin)
          allocate(sig_om(nb,nb,ndim3_omega))
          ps3=(0.d0,0.d0)
          do k=1,nqdiv_c
            v=pnt_c(:,k1)-pnt_c(:,k)
            call zone1_number(v,rb0_c,ndiv_c,k1k)
            k1k=index_k1_c(k1k)    !! for K'-K
            gtild=v-pnt_c(:,k1k)
            call ppm_factors_sig(k1,k,k1k,fif_full(1,1,1,1,ispin),
     &                           z_full(1,1,1,ispin),
     &                           a_full(1,1,1,ispin),ps3(1,1,1,k),
     &                           gtild,n_low_bnd(1,ispin))
          enddo   !! over k
          sig_om=(0.d0,0.d0)
          do ind_s2=1,ndim4_bnd(me_b+1,k1f,ispin)
            s2=n4_mpi_bnd(me_b+1,k1f,ispin)+ind_s2
            allocate(t_stat(n_pbtot_red,nbndf_bnd,nqdiv_c,2,0:n_tau/2))
            t_stat=(0.d0,0.d0)
            do ind_tau=1,ndim3_tau
              i_tau=me_t*ndim3_tau+ind_tau-1
              do it=1,2
                call tsig_stat(it,ind_tau,t_stat(1,1,1,it,i_tau),
     &                         g_tau(1,1,1,1,1,ispin),v2_q,ps3,s2,k1,
     &                         ispin,z_red(1,1,1,ispin),
     &                         a_red(1,1,1,ispin),ax(1,1,1,ispin),phf,
     &                         grsmm,grsmi,grsim,grsii,nb)
              enddo  !! over it
            enddo   !!! over ind_tau
            if(nproc_t/=1) 
     &       call dgop(t_stat,2*n_pbtot_red*nbndf_bnd*(n_tau+2)*nqdiv_c,
     &                 '  +',comm_t)
            do ind_omega=1,ndim3_omega
              i1_omega=me_t*ndim3_omega+ind_omega-1
              nom=num_omega1_omega(i1_omega)
              nomnu=num_nu_omega(i1_omega)
c              numtn=max(nomnu,n_tau+1)
              numtn=max(nom,n_tau+1)
              allocate(t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c,numtn,2))
              allocate(g_om_om(nbndf_bnd,nbndf_bnd,2,nom,npnt_c))
              g_om_om=(0.d0,0.d0)
              allocate(g_nu_om(nbndf_bnd,nbndf_bnd,2,nomnu,npnt_c))
              g_nu_om=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
              call g_interp_sig(gexa,gspl_geom,gspl_asy,gasy,g_om_om,
     &                          g_nu_om,i1_omega,nom,nomnu,ispin)
c --------------------------------------------------------------------
              if(k1==1.and.i1_omega==0.and.s2==1)
     &         call timel('**** K0_SIG starts *****************')
              call k0_sig(ind_omega,k1,k_pw,g_om_om,ps3,s2,nom,ispin,
     &                    wnu)
              if(k1==1.and.i1_omega==0.and.s2==1)
     &         call timel('**** K0_SIG finished ***************')
              do itr=2,iter_s
                t_pw=(0.d0,0.d0)
                do i_tau=0,n_tau
                  do it=1,2
                    call k_t_sig(ispin,i_tau,k1,k_pw(1,1,1,i_tau,it),
     &                           z_red(1,1,1,ispin),a_red(1,1,1,ispin),
     &                           ax(1,1,1,ispin),phf,grsmm,grsmi,grsim,
     &                           grsii,t_pw(1,1,1,i_tau+1,it))
                  enddo  !! over it
                enddo  !! over i_tau
                if(k1==1.and.i1_omega==0.and.s2==1.and.itr==2)
     &            call timel('**** K_RS and TSIG finished ********')
c --------------------------------------------------------------------
                call tsig_tau_to_omega(t_pw,ind_omega,k1,nom,ispin,
     &                                 numtn)
                if(k1==1.and.i1_omega==0.and.s2==1.and.itr==2)
     &            call timel('**** TSIG_TAU_to_OMEGA finished ****')
c --------------------------------------------------------------------
                call k_sig(ind_omega,k1,k_pw,g_om_om,ps3,nom,ispin,wnu,
     &                     numtn,t_pw,t_stat,v2_q,itr,
     &                     sig_om(1,s2,ind_omega))
                if(k1==1.and.i1_omega==0.and.s2==1.and.itr==2)
     &            call timel('**** K_SIG finished ****************')
              enddo  !! over itr
              deallocate(t_pw,g_om_om,g_nu_om)
            enddo   !! over ind_omega
            deallocate(t_stat)
            if(k1==1.and.s2==1) 
     &          call timel('**** Omega-loop finished ***********')
          enddo   !! over ind_s2
          if(nproc_b/=1) call dgop(sig_om,2*nb*nb*ndim3_omega,'  +',
     &                               comm_b)
          do ind_omega=1,ndim3_omega
            call ferm_pack_omega(sig_om(1,1,ind_omega),
     &                       sig_dyn_om(1,1,1,ind_omega,ind_k,ispin),nb,
     &                           nb,nbndf_bnd)
          enddo
          if(k1==1) call timel('**** S2-loop finished **************')
          deallocate(sig_om)
          if(k1==1) call timel('**** END of K=1 loop ***************')
        enddo  !! over ind_k
        call timel('**** END of K loop *****************')
      enddo   !! over ispin
      deallocate(v2_q,k_pw,gexa,gspl_geom,gspl_asy,gasy,g_x,ex,
     &           g_tau,t_rr_k,fif_full,z_full,
     &           a_full,wnu,g_omega,grsmm,grsmi,grsim,grsii,
     &           wnuq,a_red)
c --------------------------------------------------------------- 
      herm=0.d0
      allocate(sigt(nbndf_bnd,nbndf_bnd,2,ndim3_tau))
      do ispin=1,nspin
        do ind_k=1,ndim_k_red(me_k+1)
          k=n_mpi_k_red(me_k+1)+ind_k
          ka=k_a_from_c(k)
          n=n_low_bnd(ka,ispin)
          call omega_to_tau_baa(sig_dyn_om(1,1,1,1,ind_k,ispin),
     &                          nbndf_bnd,sigt,nbndf_bnd,n)
          do ind_tau=1,ndim3_tau
            call ferm_unpack_tau(sig_1_1(1,1,1,ind_tau,ind_k,ispin),
     &                           sigt(1,1,1,ind_tau),n,nbndf_bnd,
     &                           nbndf_bnd,3) 
            do it=1,2
              hm=check_hermit(sig_1_1(1,1,it,ind_tau,ind_k,ispin),n,
     &                        nbndf_bnd)
              herm=max(herm,hm)
            enddo
          enddo
        enddo
      enddo
      deallocate(sigt)
      if(nproc_t/=1)call dgop(herm,1,'max',comm_t)
      if(nproc_k/=1) call dgop(herm,1,'max',comm_k)
      if(maswrk) write(iun,*)' VERTEX_SIG_1 : HERM_DYN_SIGMA = ',herm
      deallocate(ps3,z_red,phf,nind_fif,ind_fif,fif0,ax)
      call timel('**** VERTEX_SIG_1_NLC finished *****')
      end
