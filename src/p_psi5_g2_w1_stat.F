      subroutine p_psi5_g2_w1_stat(tmatr,ch)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      complex*16, intent(in) :: tmatr(n_lem2_red,n_lem2_red,ndim3_nu,
     &                                ndimc_kk(me_k+1),nspin,nspin)
      integer :: k,ispin,iq,ind_nu,i_nu,i,nom,s2,i_omega,np,ind_q,nbb,
     &           ib,ind_omega,i_tau,kq,it,ir,nomnu,ind_s2,kf,
     &           ka,j,ind_k,ind_r,iq1,jspin
      real*8 :: p_herm,check_hermit,ph,de,green0,v(3),gtild(3,2),
     &          pi2
      real*8, allocatable :: poli(:,:),gexa(:,:,:,:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:,:,:),
     &                       gasy(:,:,:,:,:),t_rr(:,:,:,:,:,:,:),
     &                       g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),
     &                       fif_full(:,:,:,:,:)
      complex*16, allocatable :: k0_pw(:,:,:,:,:,:),trn(:,:),
     &                           pol(:,:,:,:),lambda_dyn(:,:,:,:,:,:),
     &                           g_om_nu(:,:,:,:,:,:),
     &                           z_red(:,:,:,:),ps3(:,:,:,:,:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           z_full(:,:,:,:),
     &                           a_full(:,:,:,:),phf(:,:),kr(:,:,:),
     &                           tn(:,:,:,:)
      pi2=pi+pi
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
      deallocate(a_red,ax,phf,ind_fif,nind_fif,fif0)
c ------ Here for zero-order vertex in full basis ------------------
      allocate(fif_full(nrel,maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ---   T-matrix in (R;tau) veriables ------------------------------
      allocate(t_rr(n_lem2_red,n_lem2_red,2,ndim3_tau,
     &              ndimc_kk(me_k+1),nspin,nspin))
      allocate(tn(n_lem2_red,n_lem2_red,ndim3_nu,ndimc_kk(me_k+1)))
      allocate(kr(nqdiv_c,n_lem2_red,n_lem2_red))
      do jspin=1,nspin
        do ispin=1,nspin
          do ind_nu=1,ndim3_nu
            kr=(0.d0,0.d0)
            do ind_q=1,ndimc_kk(me_k+1)
              iq=nc_mpi_kk(me_k+1)+ind_q
              call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,iq1)
              kr(iq1,:,:)=tmatr(:,:,ind_nu,ind_q,ispin,jspin)
            enddo
            if(nproc_k/=1) call dgop(kr,2*n_lem2_red**2*nqdiv_c,'  +',
     &                               comm_k)
            call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_lem2_red**2,kr,1)
            do ind_r=1,ndimc_kk(me_k+1)
              ir=nc_mpi_kk(me_k+1)+ind_r
              tn(:,:,ind_nu,ind_r)=kr(ir,:,:)/nqdiv_c
            enddo
          enddo   !! over ind_nu
          do ind_r=1,ndimc_kk(me_k+1)
            call nu_to_tau_c_spl(tn(1,1,1,ind_r),
     &                           t_rr(1,1,1,1,ind_r,ispin,jspin),
     &                           n_lem2_red**2)
          enddo
        enddo
      enddo
      if(nspin==1.and.irel/=2) t_rr=2.d0*t_rr
      deallocate(kr,tn)
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
      deallocate(g_tau)
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
      allocate(k0_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2,ndimc_kk(me_k+1),
     &               nspin))
c ------------------------------------------------------------------
      allocate(pol(n_pbtot_red,n_pbtot_red,0:n_nu,npnt_c))
      pol=(0.d0,0.d0)
      allocate(ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,ndimc_kk(me_k+1),
     &             nspin))
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
      do iq=1,npnt_c
        ps3=(0.d0,0.d0)
        do ispin=1,nspin
          do ind_k=1,ndimc_kk(me_k+1)
            k=nc_mpi_kk(me_k+1)+ind_k
            v=pnt_c(:,k)-pnt_c(:,iq)
            call zone1_number(v,rb0_c,ndiv_c,kq)
            kq=index_k1_c(kq)    !! for K-Q
            gtild(:,1)=0.d0
            gtild(:,2)=v-pnt_c(:,kq)
            call ppm_factors_full(k,kq,iq,fif_full(1,1,1,1,ispin),
     &                            z_full(1,1,1,ispin),
     &                            a_full(1,1,1,ispin),
     &                            ps3(1,1,1,ind_k,ispin),gtild,
     &                            n_low_bnd(1,ispin))
          enddo   !! over ind_k
        enddo   !! over ispin
        do i_nu=0,n_nu
          nom=num_omega1_adapt(i_nu)
          nomnu=num_nu2_adapt(i_nu)
          allocate(lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,
     &                        ndimc_kk(me_k+1),nspin))
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
c -------- Get zero order generalized susceptibility ---------------
            call k0_g2_t2(i_nu,iq,k0_pw,g_om_nu,ps3,s2,g_x,ex,nom)
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND finished ***************')
c ------------------------------------------------------------------
            call vertex_g2_t2(iq,i_nu,s2,nom,k0_pw,z_red,lambda_dyn,
     &                        ch,t_rr)
c ------------------------------------------------------------------
            do ispin=1,nspin
c -------------- Correction to the Polarizability ------------------
              do ind_k=1,ndimc_kk(me_k+1)
                k=nc_mpi_kk(me_k+1)+ind_k
                call dp_g2_t2(nom,i_nu,k,iq,
     &                        lambda_dyn(1,1,1,1,ind_k,ispin),
     &                        g_om_nu(1,1,1,1,1,ispin),
     &                        ps3(1,1,1,k,ispin),pol(1,s2,i_nu,iq),
     &                        ispin,trn(1,i_nu),s2)
              enddo  !! over ind_k
            enddo  !! over ispin
          enddo   !! over ind_s2
          if(iq==1.and.i_nu==0)
     &       call timel('**** S2-loop finished **************')
          deallocate(lambda_dyn,g_om_nu)
        enddo   !! over i_nu
        if(iq==1) call timel('**** END of Q=1 loop ***************')
      enddo  !! over iq
      call timel('**** END of Q loop *****************')
      deallocate(k0_pw,trn,gexa,gspl_geom,gspl_asy,gasy,g_x,ex,
     &           fif_full,z_full,a_full,t_rr)
c ---------------------------------------------------------------
      if(nproc_b/=1) call dgop(pol,2*n_pbtot_red**2*(n_nu+1)*npnt_c,
     &                         '  +',comm_b)
      if(nproc_k/=1) call dgop(pol,2*n_pbtot_red**2*(n_nu+1)*npnt_c,
     &                         '  +',comm_k)
      if(nspin==1.and.irel/=2) pol=2.d0*pol
      p_herm=0.d0
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        ka=k_a_from_c(iq)
        np=n_pbmt_red+nplwgw_red(ka)
        allocate(poli(np,np))
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          ph=check_hermit(pol(1,1,i_nu,iq),np,n_pbtot_red)
          p_herm=max(p_herm,ph)
          call pack_hermit(pol(1,1,i_nu,iq),poli,np,n_pbtot_red,np,0.d0,
     &                     1.d0)
          p_vertex_nu_q(1:np,1:np,ind_nu,ind_q)=
     &      p_vertex_nu_q(1:np,1:np,ind_nu,ind_q)+poli
        enddo
        deallocate(poli)
      enddo
      if(nproc_t/=1) call dgop(p_herm,1,'max',comm_t)
      if(nproc_k/=1) call dgop(p_herm,1,'max',comm_k)
      if(maswrk) write(iun,*)' P_PSI5_G2_W1_STAT : P_HERM = ',p_herm
      deallocate(ps3,z_red,pol)
      call timel('**** P_PSI5_G2_W1_STAT finished ****')
      end
