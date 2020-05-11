      subroutine vertex_p_vmt(iter_w,iter_t,key)
c     key = 0 - Based on DFT or HF Green's function
c     key = 1 - Based on QP Green's function
c     key = 2 - Based on general Green's function
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iter_w,iter_t,key
      logical :: spin_flips
      integer :: k,ispin,iq,ind_nu,i_nu,i,nom,it_vrt,s2,i_omega,np,nit,
     &           ind_q,nbb,ib,ind_omega,i_tau,kq,iq_line,it,nomnu,
     &           ind_s2,idum,kf,iqf,ka,j
      real*8 :: p_herm,check_hermit,ph,de,green0,v(3),gtild(3,2),
     &          lam_max,conv_1,pi2
      complex*16 :: cc
      real*8, allocatable :: poli(:,:),gexa(:,:,:,:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:,:,:),
     &                       gasy(:,:,:,:,:),g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),conv_table(:,:,:),
     &                       fif_full(:,:,:,:,:),vloc(:,:,:)
      complex*16, allocatable :: k0_pw(:,:,:,:,:,:),trn(:,:),
     &                           k_pw(:,:,:,:,:,:),pol(:,:,:,:),
     &                           lambda_dyn(:,:,:,:,:,:),
     &                           lambda_stat(:,:,:,:),
     &                           g_om_nu(:,:,:,:,:,:),
     &                           z_red(:,:,:,:),ps3(:,:,:,:,:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           vrt_dyn_old(:,:,:,:,:,:),
     &                           vrt_stat_old(:,:,:,:),
     &                           z_full(:,:,:,:),
     &                           a_full(:,:,:,:),phf(:,:)
      pi2=pi+pi
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
      deallocate(ax,phf,a_red)
c ---------------------------------------------------------------------
      call norma_wf_red
c ------ Here for zero-order vertex in full basis ------------------
      allocate(fif_full(nrel,maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ---------- Interaction matrices ----------------------------------
      allocate(vloc(n_pbmtm_red,n_pbmtm_red,natom))
      call w2_loc_stat(vloc)
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
c ------------------------------------------------------------------
      if(iter_t/=0) then
        allocate(g_rs_mm(nrel,nfun_red,nfun_red,0:n_tau,nqdiv_c,nspin))
        do ispin=1,nspin
          call g_mm_rs(g_tau(1,1,1,1,1,ispin),
     &                 g_rs_mm(1,1,1,0,1,ispin),z_red(1,1,1,ispin),
     &                 ispin)
        enddo
      endif
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
      call timel('*** VERTEX_BND_GW : G_0 finished ***')
      allocate(k0_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2,nqdiv_c,nspin))
      allocate(k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2,nqdiv_c,nspin))
      allocate(lambda_stat(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      allocate(vrt_stat_old(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      if(iter_t/=0) then
        allocate(eqr(nqdiv_c))
      endif
c ------------------------------------------------------------------
      allocate(conv_table(0:iter_w,0:7,nc_line))
      conv_table=0.d0
      allocate(pol(n_pbtot_red,n_pbtot_red,ndim3_nu,
     &             ndim_k_red(me_k+1)))
      pol=(0.d0,0.d0)
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
      call timel('**** VERTEX_P_VMT : loop starts ****')
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        if(iter_t/=0) call eqr_vertex_line(eqr,iq)
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
     &                            z_full(1,1,1,ispin),
     &                            a_full(1,1,1,ispin),
     &                            ps3(1,1,1,k,ispin),gtild,
     &                            n_low_bnd(1,ispin))
          enddo   !! over k
        enddo   !! over ispin
        lam_max=maxval(abs(ps3))
        if(iq_line/=0) conv_table(0,:,iq_line)=lam_max/nproc_t
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          nit=iter_w
          nom=num_omega1_adapt(i_nu)
          nomnu=num_nu2_adapt(i_nu)
          allocate(lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,nqdiv_c,nspin))
          allocate(vrt_dyn_old(nbndf_bnd,nbndf_bnd,nom,2,nqdiv_c,
     &                       nspin))
          if(key/=0) then
            allocate(g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,npnt_c,nspin))
            g_om_nu=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
            call g_interp_bnd(gexa,gspl_geom,gspl_asy,gasy,g_om_nu,i_nu,
     &                        nom)
          endif
c --------------------------------------------------------------------
          do ind_s2=1,ndim4_pbr(me_b+1,iq)
            s2=n4_mpi_pbr(me_b+1,iq)+ind_s2
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND starts *****************')
            lambda_dyn=(0.d0,0.d0)
c -------- Get zero order generalized susceptibility ---------------
            if(key==0) then
              call k0_bnd_exact(i_nu,iq,k0_pw,ps3,s2,g_x,ex)
            else
              call k0_bnd(ind_nu,iq,k0_pw,g_om_nu,ps3,s2,g_x,ex,nom)
            endif
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND finished ***************')
c ----------------------------------------------------------------------
            k_pw=(0.d0,0.d0)
            vrt_dyn_old=(0.d0,0.d0)
            vrt_stat_old=(0.d0,0.d0)
            do it_vrt=1,nit
              if(iq==1.and.i_nu==0.and.s2==1.and.it_vrt==1)
     &          call timel('****** Iter-n starts ***************')
              k_pw=k0_pw+k_pw
              spin_flips=.false.
              if(it_vrt<=iter_t) spin_flips=.true.
              call vertex_vmt_loop(spin_flips,iq,ind_nu,it_vrt,2,s2,nom,
     &                             nomnu,k_pw,fif0,z_red,vloc,
     &                             lambda_stat,lambda_dyn,nspin)
c --------------------------------------------------------------------
              conv_1=0.d0
              do ispin=1,nspin
                do k=1,nqdiv_c
                  do i_omega=1,nom
                    do it=1,2
                      do j=1,nbndf_bnd
                        do i=1,nbndf_bnd
                          cc=lambda_dyn(i,j,i_omega,it,k,ispin)
     &                      +lambda_stat(i,j,k,ispin)
     &                      -vrt_dyn_old(i,j,i_omega,it,k,ispin)
     &                      -vrt_stat_old(i,j,k,ispin)
                          conv_1=max(conv_1,abs(cc))
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              if(i_nu<=7.and.iq_line/=0) then
                conv_table(it_vrt,i_nu,iq_line)=
     &              max(conv_1,conv_table(it_vrt,i_nu,iq_line))
              endif
              vrt_dyn_old=lambda_dyn
              vrt_stat_old=lambda_stat
              do ispin=1,nspin
                if(it_vrt<nit) then
                  k_pw(:,:,:,:,:,ispin)=(0.d0,0.d0)
                  do k=1,nqdiv_c
c -------- Get correction to the generalized susceptibility --------
                    if(key==0) then
                      call k_bnd_exact(ispin,nom,ind_nu,k,iq,
     &                                 lambda_dyn(1,1,1,1,k,ispin),
     &                                 k_pw(1,1,0,1,k,ispin),
     &                                 lambda_stat(1,1,k,ispin),
     &                                 g_x(1,0,1,ispin),ex(1,1,ispin))
                    else
                      call k_bnd(ispin,nom,ind_nu,k,iq,
     &                           lambda_dyn(1,1,1,1,k,ispin),
     &                           k_pw(1,1,0,1,k,ispin),
     &                           lambda_stat(1,1,k,ispin),
     &                           g_om_nu(1,1,1,1,1,ispin),
     &                           g_x(1,0,1,ispin),ex(1,1,ispin))
                    endif
                  enddo  !! over k
                else if(it_vrt==nit) then
c -------------- Correction to the Polarizability ------------------
                  do k=1,nqdiv_c
                    if(key==0) then
                      call dp_bnd_exact(nom,i_nu,k,iq,
     &                                  lambda_dyn(1,1,1,1,k,ispin),
     &                                  lambda_stat(1,1,k,ispin),
     &                                  ps3(1,1,1,k,ispin),
     &                                  pol(1,s2,ind_nu,ind_q),ispin,
     &                                  g_x(1,0,1,ispin),ex(1,1,ispin),
     &                                  trn(1,i_nu),1,s2)
                    else
                      call dp_bnd(nom,ind_nu,k,iq,
     &                            lambda_dyn(1,1,1,1,k,ispin),
     &                            lambda_stat(1,1,k,ispin),
     &                            g_om_nu(1,1,1,1,1,ispin),
     &                            ps3(1,1,1,k,ispin),
     &                            pol(1,s2,ind_nu,ind_q),ispin,
     &                            g_x(1,0,1,ispin),ex(1,1,ispin),
     &                            trn(1,i_nu),1,s2)
                    endif
                  enddo  !! over k
                endif
              enddo  !! over ispin
            enddo  !! over it_vrt
          enddo   !! over ind_s2
          if(iq==1.and.i_nu==0)
     &       call timel('**** S2-loop finished **************')
          deallocate(lambda_dyn,vrt_dyn_old)
          if(key/=0) deallocate(g_om_nu)
        enddo   !! over ind_nu
        if(iq==1) call timel('**** Nu-loop finished **************')
        if(iq==1) call timel('**** END of Q=1 loop ***************')
      enddo  !! over ind_q
      call timel('**** END of Q loop *****************')
      deallocate(k0_pw,k_pw,trn,lambda_stat,gexa,gspl_geom,gspl_asy,
     &           gasy,g_x,ex,vrt_stat_old,fif_full,z_full,a_full,vloc)
      if(iter_t/=0) then
        deallocate(g_rs_mm,eqr)
      endif
c ---------------------------------------------------------------
      if(nproc_b/=1) call dgop(pol,
     &                2*n_pbtot_red**2*ndim3_nu*ndim_k_red(me_k+1),
     &                         '  +',comm_b)
      if(nspin==1.and.irel/=2) pol=2.d0*pol
      p_herm=0.d0
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        ka=k_a_from_c(iq)
        np=n_pbmt_red+nplwgw_red(ka)
        allocate(poli(np,np))
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          ph=check_hermit(pol(1,1,ind_nu,ind_q),np,n_pbtot_red)
          p_herm=max(p_herm,ph)
          call pack_hermit(pol(1,1,ind_nu,ind_q),poli,np,
     &                     n_pbtot_red,np,0.d0,1.d0)
          p_vertex_nu_q(1:np,1:np,ind_nu,ind_q)=
     &      p_vertex_nu_q(1:np,1:np,ind_nu,ind_q)+poli
        enddo
        deallocate(poli)
      enddo
      if(nproc_t/=1) call dgop(p_herm,1,'max',comm_t)
      if(nproc_k/=1) call dgop(p_herm,1,'max',comm_k)
      if(maswrk) write(iun,*)' VERTEX_P_VMT : P_HERM = ',p_herm
      call output_p_bnd_nu(1,pol,n_pbtot_red)
      deallocate(ps3,z_red,nind_fif,ind_fif,fif0)
c ------------------------------------------------------------------
      idum=(iter_w+1)*8*nc_line
      if(nproc_b/=1) call dgop(conv_table,idum,'max',comm_b)
      if(nproc_k/=1) call dgop(conv_table,idum,'  +',comm_k)
      if(nproc_t/=1) call dgop(conv_table,idum,'  +',comm_t)
      call output_bs_conv(conv_table,iter_w)
      deallocate(pol,conv_table)
      call timel('**** VERTEX_P_VMT finished *********')
      end
