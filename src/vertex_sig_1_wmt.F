      subroutine vertex_sig_1_wmt(sig_1_1)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(out) :: sig_1_1(nbndf_bnd,nbndf_bnd,2,
     &                                   ndim3_tau,ndim_k_red(me_k+1),
     &                                   nspin)
      integer :: k,ispin,i_nu,i,nom,s2,i_omega,
     &           ind_omega,mx2,i_tau,it,ind_tau,n,nomnu,
     &           ind_k,ka,j,i1_omega,ind_nu,
     &           ir,numtn,i1_tau,ind_ir0,iatom,ir0,indi,indj,irl,
     &           jsort,jatom,isort,n1,np,np1
      real*8 :: w1,pib,pi2,phase
      complex*16 :: cc
      real*8, allocatable :: grsmm(:,:,:,:,:),
     &                       wexa(:,:,:),fifim(:,:,:,:,:),
     &                       wspl_geom(:,:,:,:),wspl_asy(:,:,:,:),
     &                       wasy(:,:),g0(:,:),
     &                       g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),
     &                       wnuq(:,:,:,:),wloc(:,:,:,:),
     &                       vloc(:,:,:),gexa(:,:,:,:,:),
     &                       gspl_geom(:,:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:,:),
     &                       gasy(:,:,:,:),fif(:,:,:,:,:),
     &                       sig_dyn(:,:,:,:,:,:),
     &                       sigk(:,:,:,:,:)
      complex*16, allocatable :: k_pw(:,:,:,:),gtau(:,:,:),
     &                           t_pw(:,:,:,:),gom(:,:,:),
     &                           g_om_om(:,:,:,:),g_nu_om(:,:,:,:),
     &                           z_red(:,:,:,:),sigt(:,:,:,:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           sig_stat(:,:,:,:,:),
     &                           sig_semi(:,:,:,:,:),
     &                           sig_om(:,:),phf(:,:),
     &                           sig_dyn_om(:,:,:,:),sigr(:,:,:),
     &                           tmp(:,:)
      pi2=pi+pi
      pib=pi/betta_t
      mx2=maxel_red**2
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
      deallocate(a_red,phf,ax,nind_fif,ind_fif,fif0)
c ---------------------------------------------------------------
      allocate(fifim(nrel,maxel_red,maxel_red,n_pbmtm_red,nsort))
      allocate(fif(nrel,maxel_red,maxel_red,n_pbmtm_red,nsort))
c ---------- Interaction matrices in product basis ----
      allocate(wloc(n_pbmtm_red,n_pbmtm_red,0:n_nu,natom))
      allocate(vloc(n_pbmtm_red,n_pbmtm_red,natom))
      call w_loc_from_red(vloc,wloc)
c ------ W-interpolation - on regullar NU for Semi-Sigma ---------
      allocate(wnuq(n_pbmtm_red,n_pbmtm_red,ndim3_nu,natom))
      allocate(wexa(n_pbmtm_red,n_pbmtm_red,0:n_nu_exa))
      allocate(wspl_geom(n_pbmtm_red,n_pbmtm_red,0:n_nu_geom+2,4))
      allocate(wspl_asy(n_pbmtm_red,n_pbmtm_red,n_nu_asy+1,4))
      allocate(wasy(n_pbmtm_red,n_pbmtm_red))
      allocate(tmp(n_pbmtm_red,n_pbmtm_red))
      do iatom=1,natom
        isort=is(iatom)
        call w_spl_loc_red(iatom,wloc(1,1,0,iatom),wexa,wspl_geom,
     &                     wspl_asy,wasy)
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          w1=w_nu(i_nu)
          call interp_w_nu_re_spl(wexa,wspl_geom,wspl_asy,wasy,
     &                            tmp,w1,n_pbmt0_red(isort),n_pbmtm_red)
          wnuq(:,:,ind_nu,iatom)=tmp
        enddo
      enddo
      deallocate(wexa,wspl_geom,wspl_asy,wasy,tmp)
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
      deallocate(g_omega)
c -----------------
      allocate(grsmm(nrel,nfun_red,nfun_red,0:n_tau,nqdiv_c))
c -------------	
      allocate(gexa(nfun_red,nfun_red,2,0:n_omega_exa,nqdiv_c))
      allocate(gspl_geom(nfun_red,nfun_red,2,0:n_omega_geom+2,4,
     &                   nqdiv_c))
      allocate(gspl_asy(nfun_red,nfun_red,2,n_omega_asy+1,4,nqdiv_c))
      allocate(gasy(nfun_red,nfun_red,2,nqdiv_c))
      call timel('*** VERTEX_SIG_1 : G_0 finished ****')
c -------------------------------------------------------------------
      allocate(sig_dyn_om(maxel_red,maxel_red,nrr_max_c,ndim3_omega))
      allocate(sig_dyn(nfun_red,nfun_red,2,ndim3_tau,npnt_c,nspin))
      allocate(sig_semi(maxel_red,maxel_red,nrr_max_c,2,ndim3_tau))
      allocate(sig_stat(maxel_red,maxel_red,nrr_max_c,2,ndim3_tau))
      call timel('*** most allocations finished ******')
c ------------------------------------------------------------------
      do ispin=1,nspin
        call fi0_full_red(fifim,ispin,ispin)
        do isort=1,nsort
          do j=1,lfunm_red(isort)
            do i=1,lfunm_red(isort)
              do irl=1,nrel
                fif(irl,i,j,:,isort)=fifim(irl,i,j,:,isort)
              enddo
            enddo
          enddo
        enddo
c  -----------------------------------------------------------------
        call g_mm_rs(g_tau(1,1,1,1,1,ispin),grsmm,z_red(1,1,1,ispin),
     &               ispin)
        call timel('**** G_REAL_SPACE finished *********')
        gexa=0.d0
        gspl_geom=0.d0
        gspl_asy=0.d0
        gasy=0.d0
c --------- For Interpolation of G(w') ------------------------------
        allocate(gtau(nfun_red,nfun_red,0:n_tau))
        allocate(gom(nfun_red,nfun_red,0:n_omega))
        allocate(g0(0:n_omega,2))
        do ir=1,nqdiv_c
          do j=1,nfun_red
            do i=1,nfun_red
              if(irel/=2) then
                gtau(i,j,:)=grsmm(1,i,j,:,ir)
              else
                gtau(i,j,:)=dcmplx(grsmm(1,i,j,:,ir),
     &                             grsmm(2,i,j,:,ir))
              endif
            enddo
          enddo
          call tau_to_omega_c(gom,gtau,nfun_red**2)
          do j=1,nfun_red
            do i=1,nfun_red
              g0=0.d0
              do i_omega=0,n_omega
                g0(i_omega,1)=real(gom(i,j,i_omega))
                g0(i_omega,2)=imag(gom(i,j,i_omega))
              enddo
              do it=1,2
                call spline_inhmg(w_omega(n_omega_exa-1),
     &                            g0(n_omega_exa-1:,it),
     &                            gspl_geom(i,j,it,:,1,ir),
     &                            gspl_geom(i,j,it,:,2,ir),
     &                            gspl_geom(i,j,it,:,3,ir),
     &                            gspl_geom(i,j,it,:,4,ir),
     &                            n_omega_geom+3,0,0.d0,0.d0)
                call spline_inhmg(xm_omega(0),
     &                            g0(n_omega_exa+n_omega_geom:,it),
     &                            gspl_asy(i,j,it,:,1,ir),
     &                            gspl_asy(i,j,it,:,2,ir),
     &                            gspl_asy(i,j,it,:,3,ir),
     &                            gspl_asy(i,j,it,:,4,ir),
     &                            n_omega_asy+1,0,0.d0,0.d0)
                gasy(i,j,it,ir)=g0(n_omega,it)
                gexa(i,j,it,:,ir)=g0(0:n_omega_exa,it)
              enddo
            enddo
          enddo
        enddo
        deallocate(g0,gtau,gom)
        call timel('**** VERTEX_SIG_1 : loop starts ****')
        sig_dyn_om=(0.d0,0.d0)
        sig_stat=(0.d0,0.d0)
        sig_semi=(0.d0,0.d0)
        do ind_ir0=1,ndim_nrr_c(me_k+1)
          ir0=n_mpi_nrr_c(me_k+1)+ind_ir0
          ir=nrr_red_c(1,ir0)
          iatom=nrr_red_c(2,ir0)
          jatom=nrr_red_c(3,ir0)
          jsort=is(jatom)
          isort=is(iatom)
          n=lfunm_red(isort)
          n1=lfunm_red(jsort)
          np=n_pbmt0_red(isort)
          np1=n_pbmt0_red(jsort)
          indi=io_lem_red(iatom)
          indj=io_lem_red(jatom)
          allocate(sig_om(n,n1))
          allocate(k_pw(n,np1,0:n_tau,2))
          do ind_omega=1,ndim3_omega
            i1_omega=me_t*ndim3_omega+ind_omega-1
            nom=num_omega1_omega(i1_omega)
            nomnu=num_nu_omega(i1_omega)
            numtn=max(nomnu,n_tau+1)
            allocate(t_pw(np,n1,numtn,2))
            allocate(g_om_om(n,n1,2,nom))
            g_om_om=(0.d0,0.d0)
            allocate(g_nu_om(n,n1,2,nomnu))
            g_nu_om=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
            call g_interp_wmt(gexa,gspl_geom,gspl_asy,gasy,g_om_om,
     &                        g_nu_om,i1_omega,nom,nomnu,ispin,n,n1,ir,
     &                        indi,indj)
c --------------------------------------------------------------------
            sig_om=(0.d0,0.d0)
            do s2=1,n1
              if(ir0==1.and.i1_omega==0.and.s2==1)
     &         call timel('**** K_SIG starts ******************')
              call k0_sig_wmt(ind_omega,jatom,k_pw,g_om_om,n,n1,
     &                        np1,fif(1,1,1,1,jsort),s2,nom,ispin,
     &                        wloc(1,1,0,jatom))
              if(ir0==1.and.i1_omega==0.and.s2==1)
     &         call timel('**** K_SIG finished ****************')
              t_pw=(0.d0,0.d0)
              do i_tau=0,n_tau
                i1_tau=n_tau-i_tau
                do it=1,2
                  if(irel/=2) then
                    call tsig_wmt(grsmm(1,1,1,0,ir),k_pw(1,1,i_tau,it),
     &                            t_pw(1,1,i_tau+1,it),fif,i1_tau,n,n1,
     &                            np,np1,iatom,jatom)
                  else
                    call tsig_wmt_c(grsmm(1,1,1,0,ir),
     &                              k_pw(1,1,i_tau,it),
     &                              t_pw(1,1,i_tau+1,it),fif,i1_tau,n,
     &                              n1,np,np1,iatom,jatom)
                  endif
                enddo  !! over it
              enddo  !! over i_tau
              if(ir0==1.and.i1_omega==0.and.s2==1)
     &            call timel('**** TSIG_WMT finished *************')
c --------------------------------------------------------------------
              call tsig_tau_to_nu_wmt(t_pw,ind_omega,nomnu,ispin,
     &                                numtn,np,n1)
              if(ir0==1.and.i1_omega==0.and.s2==1)
     &          call timel('**** TSIG_TAU_to_NU finished *******')
c --------------------------------------------------------------------
              allocate(wexa(n_pbmtm_red,n_pbmtm_red,0:n_nu_exa))
              allocate(wspl_geom(n_pbmtm_red,n_pbmtm_red,0:n_nu_geom+2,
     &                           4))
              allocate(wspl_asy(n_pbmtm_red,n_pbmtm_red,n_nu_asy+1,4))
              allocate(wasy(n_pbmtm_red,n_pbmtm_red))
              call w_spl_loc_red(iatom,wloc(1,1,0,iatom),wexa,wspl_geom,
     &                           wspl_asy,wasy)
              call tsig_end_dyn_wmt(nomnu,ind_omega,t_pw,g_nu_om,
     &                              sig_om(1,s2),wexa,wspl_geom,
     &                              wspl_asy,wasy,numtn,fif,n,n1,np,
     &                              iatom)
              deallocate(wexa,wspl_geom,wspl_asy,wasy)
              if(ir0==1.and.i1_omega==0.and.s2==1)
     &          call timel('**** TSIG_END_DYN_WMT finished *****')
            enddo   !! over s2
            if(ir0==1.and.i1_omega==0)
     &       call timel('**** S2-loop finished **************')
            deallocate(t_pw,g_om_om,g_nu_om)
            sig_dyn_om(1:n,1:n1,ir0,ind_omega)=sig_om
          enddo   !! over ind_omega
          deallocate(k_pw,sig_om)
          if(ir0==1) call timel('**** Omega-loop finished ***********')
c ------- Static and Semidynamic parts -------------------------------
          allocate(t_pw(np,n1,2,ndim3_tau))
          do s2=1,n1
            t_pw=(0.d0,0.d0)
            do ind_tau=1,ndim3_tau
              do it=1,2
                call sig_1_stat_wmt(it,ind_tau,t_pw(1,1,it,ind_tau),
     &                             vloc,fif,s2,ispin,grsmm(1,1,1,0,ir),
     &                             sig_stat(1,s2,ir0,it,ind_tau),n,n1,
     &                             np,np1,iatom,jatom,ir)
              enddo  !! over it
            enddo   !!! over ind_tau
            call sig_1_semi_wmt(t_pw,grsmm(1,1,1,0,ir),
     &                          wnuq(1,1,1,iatom),fif(1,1,1,1,isort),s2,
     &                          ispin,sig_semi,ir0,n,n1,np,indi,indj)
          enddo   !! over s2
          deallocate(t_pw)
          if(ir0==1) call timel('**** END of IR0=1 loop *************')
        enddo  !! over ind_ir0
        if(nproc_k/=1) then
          call dgop(sig_dyn_om,2*maxel_red**2*ndim3_omega*nrr_max_c,
     &              '  +',comm_k)
          call dgop(sig_stat,4*maxel_red**2*ndim3_tau*nrr_max_c,
     &              '  +',comm_k)
          call dgop(sig_semi,4*maxel_red**2*ndim3_tau*nrr_max_c,
     &              '  +',comm_k)
        endif
        allocate(sigk(nfun_red,nfun_red,2,ndim3_omega,nqdiv_c))
        allocate(sigr(nfun_red,nfun_red,nqdiv_c))
        allocate(sig_om(nfun_red,nfun_red))
        do ind_omega=1,ndim3_omega
          if(irel/=2) call sym_sig_rs_c(sig_dyn_om(1,1,1,ind_omega),
     &                                  sigr)
          if(irel==2) call sym_sig_rs_c_rel(sig_dyn_om(1,1,1,ind_omega),
     &                                      sigr)
          do k=1,npnt_c
            sig_om=(0.d0,0.d0)
            do ir=1,nqdiv_c
              phase=pi2*dot_product(pnt_c(:,k),rr_coarse(:,ir))
              cc=dcmplx(cos(phase),sin(phase))
              sig_om=sig_om+cc*sigr(:,:,ir)
            enddo
            call ferm_pack_omega(sig_om,sigk(1,1,1,ind_omega,k),
     &                           nfun_red,nfun_red,nfun_red)
          enddo
        enddo
        deallocate(sigr,sig_om)
        do k=1,npnt_c
          call omega_to_tau_baa(sigk(1,1,1,1,k),nfun_red,
     &                          sig_dyn(1,1,1,1,k,ispin),nfun_red,
     &                          nfun_red)
        enddo
        deallocate(sigk)
        allocate(sigr(nfun_red,nfun_red,nqdiv_c))
        allocate(sigt(nfun_red,nfun_red,2,npnt))
        do ind_tau=1,ndim3_tau
          sigt=(0.d0,0.d0)
          do it=1,2
            if(irel/=2) call sym_sig_rs_c(sig_stat(1,1,1,it,ind_tau),
     &                                    sigr)
            if(irel==2) call sym_sig_rs_c_rel(sig_stat(1,1,1,it,
     &                                                 ind_tau),sigr)
            do k=1,npnt_c
              do ir=1,nqdiv_c
                phase=pi2*dot_product(pnt_c(:,k),rr_coarse(:,ir))
                cc=dcmplx(cos(phase),sin(phase))
                sigt(:,:,it,k)=sigt(:,:,it,k)+cc*sigr(:,:,ir)
              enddo
            enddo
          enddo   !! over it
          do k=1,npnt_c
            call ferm_pack_tau(sigt(1,1,1,k),
     &                         sig_dyn(1,1,1,ind_tau,k,ispin),
     &                         nfun_red,nfun_red,nfun_red)
          enddo
          sigt=(0.d0,0.d0)
          do it=1,2
            if(irel/=2) call sym_sig_rs_c(sig_semi(1,1,1,it,ind_tau),
     &                                    sigr)
            if(irel==2) call sym_sig_rs_c_rel(sig_semi(1,1,1,it,
     &                                                 ind_tau),sigr)
            do k=1,npnt_c
              do ir=1,nqdiv_c
                phase=pi2*dot_product(pnt_c(:,k),rr_coarse(:,ir))
                cc=dcmplx(cos(phase),sin(phase))
                sigt(:,:,it,k)=sigt(:,:,it,k)+cc*sigr(:,:,ir)
              enddo
            enddo
          enddo   !! over it
          do k=1,npnt_c
            do it=1,2
              sigt(:,:,it,k)=sigt(:,:,it,k)
     &                      +conjg(transpose(sigt(:,:,it,k)))
            enddo
            call ferm_pack_tau(sigt(1,1,1,k),
     &                         sig_dyn(1,1,1,ind_tau,k,ispin),
     &                         nfun_red,nfun_red,nfun_red)
          enddo
        enddo
        deallocate(sigr,sigt)
        call timel('**** END of IR0 loop ***************')
      enddo   !! over ispin
      deallocate(gexa,gspl_geom,gspl_asy,gasy,g_tau,fifim,fif,
     &           grsmm,wnuq,wloc,vloc,sig_dyn_om,sig_semi,sig_stat)
c ---------------------------------------------------------------
      allocate(sig_om(nfun_red,nfun_red))
      do ispin=1,nspin
        do ind_k=1,ndim_k_red(me_k+1)
          k=n_mpi_k_red(me_k+1)+ind_k
          ka=k_a_from_c(k)
          n=n_low_bnd(ka,ispin)
          allocate(tmp(n,nfun_red))
          do ind_tau=1,ndim3_tau
            do it=1,2
              call ferm_unpack_tau(sig_om,
     &                             sig_dyn(1,1,1,ind_tau,k,ispin),
     &                             nfun_red,nfun_red,nfun_red,it)
              call zgemm('c','n',n,nfun_red,nfun_red,(1.d0,0.d0),
     &                   z_red(1,1,k,ispin),nfun_red,sig_om,nfun_red,
     &                   (0.d0,0.d0),tmp,n)
              call zgemm('n','n',n,n,nfun_red,(1.d0,0.d0),tmp,n,
     &                   z_red(1,1,k,ispin),nfun_red,(0.d0,0.d0),
     &                   sig_1_1(1,1,it,ind_tau,ind_k,ispin),nbndf_bnd)
            enddo
          enddo
          deallocate(tmp)
        enddo
      enddo
      deallocate(sig_dyn,z_red,sig_om)
      call timel('**** VERTEX_SIG_1_WMT finished *****')
      end
