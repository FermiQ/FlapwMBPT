      subroutine vertex_sig_1_loc(sig_1_1)
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
      integer :: k,ispin,i_nu,i,nom,s2,i_omega,nb,nbb,
     &           ind_omega,mx2,i_tau,it,ind_tau,n,nomnu,
     &           ind_k,ka,j,i1_omega,ind_nu,
     &           numtn,i1_tau,iatom,indi,irl,isort,np
      real*8 :: w1,pib,pi2
      real*8, allocatable :: v00(:,:),w00(:,:,:),
     &                       wexa(:,:,:),fifim(:,:,:,:,:),
     &                       wspl_geom(:,:,:,:),wspl_asy(:,:,:,:),
     &                       wasy(:,:),g0(:,:),wsem(:,:,:),
     &                       g_tau(:,:,:,:,:,:,:),
     &                       wnuq(:,:,:,:),wloc(:,:,:,:),
     &                       vloc(:,:,:),gexa(:,:,:,:),
     &                       gspl_geom(:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:),
     &                       gasy(:,:,:),fif(:,:,:,:),
     &                       sig_dyn(:,:,:,:),gtau(:,:,:,:),
     &                       sig_dyn_om(:,:,:,:)
      complex*16, allocatable :: k_pw(:,:,:,:),g_omega(:,:,:,:,:),
     &                           t_pw(:,:,:,:),ax(:,:,:,:),
     &                           g_om_om(:,:,:,:),g_nu_om(:,:,:,:),
     &                           z_red(:,:,:,:),
     &                           a_red(:,:,:,:),
     &                           sig_stat(:,:,:,:),sig_semi(:,:,:,:),
     &                           sig_om(:,:),phf(:,:),
     &                           tmp(:,:),tmp1(:,:)
      pi2=pi+pi
      pib=pi/betta_t
      mx2=maxel_red**2
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
      deallocate(a_red,ax,phf,nind_fif,ind_fif,fif0)
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
      allocate(g_tau(nrel,maxel_red,maxel_red,2,ndim3_tau,nspin,natom))
      allocate(g_omega(maxel_red,maxel_red,ndim3_omega,nspin,natom))
      call g_loc_from_full(g_omega,g_tau,z_red)
      allocate(fifim(nrel,maxel_red,maxel_red,n_pbmtm_red,nsort))
c ------------------------------------------------------------------
      sig_1_1=(0.d0,0.d0)
      do ispin=1,nspin
        call fi0_full_red(fifim,ispin,ispin)
        do iatom=1,natom
          isort=is(iatom)
          n=lfunm_red(isort)
          if(n==0) cycle
          np=n_pbmt0_red(isort)
          indi=io_lem_red(iatom)
          allocate(fif(nrel,n,n,np))
c ------------------------------------------------------------------
          do j=1,n
            do i=1,n
              do irl=1,nrel
                fif(irl,i,j,1:np)=fifim(irl,i,j,1:np,isort)
              enddo
            enddo
          enddo
c  -----------------------------------------------------------------
          allocate(w00(np,np,0:n_nu))
          allocate(v00(np,np))
          v00=vloc(1:np,1:np,iatom)
          w00=wloc(1:np,1:np,:,iatom)
          allocate(wsem(np,np,ndim3_nu))
          wsem=wnuq(1:np,1:np,:,iatom)
c  -----------------------------------------------------------------
          allocate(gtau(nrel,n,n,0:n_tau))
          gtau=0.d0
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=me_t*ndim3_tau+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
              gtau(:,:,:,i_tau)=
     &          g_tau(:,1:n,1:n,it,ind_tau,ispin,iatom)
            enddo
          enddo
          if(nproc_t/=1) call dgop(gtau,n*n*nrel*(n_tau+1),'  +',
     &                               comm_t)
c  -----------------------------------------------------------------
          allocate(gexa(n,n,2,0:n_omega_exa))
          allocate(gspl_geom(n,n,2,0:n_omega_geom+2,4))
          allocate(gspl_asy(n,n,2,n_omega_asy+1,4))
          allocate(gasy(n,n,2))
          gexa=0.d0
          gspl_geom=0.d0
          gspl_asy=0.d0
          gasy=0.d0
c --------- For Interpolation of G(w') ---------------------------------
          nbb=n**2
          allocate(g0(0:n_omega,2))
          do j=1,n
            do i=1,n
              g0=0.d0
              do ind_omega=1,ndim3_omega
                i_omega=me_t*ndim3_omega+ind_omega-1
                g0(i_omega,1)=real(g_omega(i,j,ind_omega,ispin,iatom))
                g0(i_omega,2)=imag(g_omega(i,j,ind_omega,ispin,iatom))
              enddo
              if(nproc_t/=1) call dgop(g0,2*(n_omega+1),'  +',
     &                                   comm_t)
              do it=1,2
                call spline_inhmg(w_omega(n_omega_exa-1),
     &                            g0(n_omega_exa-1:,it),
     &                            gspl_geom(i,j,it,:,1),
     &                            gspl_geom(i,j,it,:,2),
     &                            gspl_geom(i,j,it,:,3),
     &                            gspl_geom(i,j,it,:,4),
     &                            n_omega_geom+3,0,0.d0,0.d0)
                call spline_inhmg(xm_omega(0),
     &                            g0(n_omega_exa+n_omega_geom:,it),
     &                            gspl_asy(i,j,it,:,1),
     &                            gspl_asy(i,j,it,:,2),
     &                            gspl_asy(i,j,it,:,3),
     &                            gspl_asy(i,j,it,:,4),
     &                            n_omega_asy+1,0,0.d0,0.d0)
                gasy(i,j,it)=g0(n_omega,it)
                gexa(i,j,it,:)=g0(0:n_omega_exa,it)
              enddo
            enddo
          enddo
          deallocate(g0)
c  -----------------------------------------------------------------
          allocate(sig_dyn_om(n,n,2,ndim3_omega))
          allocate(sig_dyn(n,n,2,ndim3_tau))
          allocate(sig_semi(n,n,2,ndim3_tau))
          allocate(sig_stat(n,n,2,ndim3_tau))
          sig_dyn_om=0.d0
          sig_stat=(0.d0,0.d0)
          sig_semi=(0.d0,0.d0)
          allocate(sig_om(n,n))
          allocate(k_pw(n,np,0:n_tau,2))
          do ind_omega=1,ndim3_omega
            i1_omega=me_t*ndim3_omega+ind_omega-1
            nom=num_omega1_omega(i1_omega)
            nomnu=num_nu_omega(i1_omega)
            numtn=max(nomnu,n_tau+1)
            allocate(t_pw(np,n,numtn,2))
            allocate(g_om_om(n,n,2,nom))
            g_om_om=(0.d0,0.d0)
            allocate(g_nu_om(n,n,2,nomnu))
            g_nu_om=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
            call g_interp_loc_sig(gexa,gspl_geom,gspl_asy,gasy,g_om_om,
     &                            g_nu_om,i1_omega,nom,nomnu,n)
c --------------------------------------------------------------------
            sig_om=(0.d0,0.d0)
            do s2=1,n
              call k0_sig_loc(ind_omega,k_pw,g_om_om,n,np,fif,s2,nom,
     &                        w00)
              t_pw=(0.d0,0.d0)
              do i_tau=0,n_tau
                i1_tau=n_tau-i_tau
                do it=1,2
                  if(irel/=2) then
                    call tsig_loc(gtau(1,1,1,i1_tau),k_pw(1,1,i_tau,it),
     &                            t_pw(1,1,i_tau+1,it),fif,n,np)
                  else
                    call tsig_loc_c(gtau(1,1,1,i1_tau),
     &                              k_pw(1,1,i_tau,it),
     &                              t_pw(1,1,i_tau+1,it),fif,n,np)
                  endif
                enddo  !! over it
              enddo  !! over i_tau
c --------------------------------------------------------------------
              call tsig_tau_to_nu_wmt(t_pw,ind_omega,nomnu,ispin,numtn,
     &                                np,n)
c --------------------------------------------------------------------
              allocate(wexa(np,np,0:n_nu_exa))
              allocate(wspl_geom(np,np,0:n_nu_geom+2,4))
              allocate(wspl_asy(np,np,n_nu_asy+1,4))
              allocate(wasy(np,np))
              call w_spl_loc(np,np,np,w00,wexa,wspl_geom,wspl_asy,wasy)
              call tsig_end_dyn_loc(nomnu,ind_omega,t_pw,g_nu_om,
     &                              sig_om(1,s2),wexa,wspl_geom,
     &                              wspl_asy,wasy,numtn,fif,n,np)
              deallocate(wexa,wspl_geom,wspl_asy,wasy)
            enddo   !! over s2
            deallocate(t_pw,g_om_om,g_nu_om)
            call ferm_pack_omega(sig_om,sig_dyn_om(1,1,1,ind_omega),n,n,
     &                           n)
          enddo   !! over ind_omega
          deallocate(k_pw,sig_om,gexa,gspl_geom,gspl_asy,gasy,w00)
c ------- Static and Semidynamic parts -------------------------------
          allocate(t_pw(np,n,2,ndim3_tau))
          do s2=1,n
            t_pw=(0.d0,0.d0)
            do ind_tau=1,ndim3_tau
              do it=1,2
                call sig_1_stat_loc(it,ind_tau,t_pw(1,1,it,ind_tau),
     &                              v00,fif,s2,gtau,
     &                              sig_stat(1,s2,it,ind_tau),n,np)
              enddo  !! over it
            enddo   !!! over ind_tau
            call sig_1_semi_loc(t_pw,gtau,wsem,fif,s2,sig_semi,n,np)
          enddo   !! over s2
          deallocate(t_pw,gtau,v00,fif)
          call omega_to_tau_baa(sig_dyn_om,n,sig_dyn,n,n)
          deallocate(sig_dyn_om)
          allocate(tmp(n,n))
          allocate(tmp1(nbndf_bnd,n))
          do ind_tau=1,ndim3_tau
            do it=1,2
              call ferm_unpack_tau(tmp,sig_dyn(1,1,1,ind_tau),n,n,n,it)
              tmp=tmp+sig_stat(:,:,it,ind_tau)+sig_semi(:,:,it,ind_tau)
     &                +conjg(transpose(sig_semi(:,:,it,ind_tau)))
              do ind_k=1,ndim_k_red(me_k+1)
                k=n_mpi_k_red(me_k+1)+ind_k
                ka=k_a_from_c(k)
                nb=n_low_bnd(ka,ispin)
                call zgemm('c','n',nb,n,n,(1.d0,0.d0),
     &                     z_red(indi,1,k,ispin),nfun_red,tmp,n,
     &                     (0.d0,0.d0),tmp1,nbndf_bnd)
                call zgemm('n','n',nb,nb,n,(1.d0,0.d0),tmp1,nbndf_bnd,
     &                     z_red(indi,1,k,ispin),nfun_red,(1.d0,0.d0),
     &                     sig_1_1(1,1,it,ind_tau,ind_k,ispin),
     &                     nbndf_bnd)
              enddo
            enddo
          enddo
          deallocate(tmp,tmp1,sig_dyn,sig_semi,sig_stat,wsem)
        enddo  !! over iatom
      enddo   !! over ispin
      deallocate(fifim,wnuq,wloc,vloc,g_tau,g_omega,z_red)
      call timel('**** VERTEX_SIG_1_LOC finished *****')
      end
