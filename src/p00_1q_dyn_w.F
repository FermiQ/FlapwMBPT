      subroutine p00_1q_dyn_w(pw,ew,e1w,rw)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(out) :: pw(0:n_nu),ew(0:n_nu),e1w(0:n_nu),
     &                           rw(0:n_nu)
      logical :: spin_flips
      integer :: k,ispin,iq,ind_nu,i_nu,i,nom,it_vrt,s2,i_omega,np,nit,
     &           nbb,ib,ind_omega,mx2,i_tau,kq,kst,it,
     &           nomnu,ind_s2,idum,kf,iqf,ka,j,nrr3
      real*8 :: p_herm,check_hermit,ph,de,green0,v(3),gtild(3,2),
     &          lam_max,conv_1
      complex*16 :: cc
      integer, allocatable :: indrr3(:,:,:)
      real*8, allocatable :: v4_r(:,:,:),w4_r(:,:,:,:),
     &                       gexa(:,:,:,:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:,:,:),
     &                       gasy(:,:,:,:,:),
     &                       w_tau(:,:,:,:),g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),vmi(:,:,:),
     &                       wmi(:,:,:,:),vim(:,:,:),
     &                       wim(:,:,:,:),vii(:),
     &                       wii(:,:),conv_table(:,:),
     &                       fif_full(:,:,:,:)
      complex*16, allocatable :: k0_pw(:,:,:,:,:,:),trn(:,:),
     &                           k_pw(:,:,:,:,:,:),pol(:,:,:),
     &                           lambda_dyn(:,:,:,:,:,:),
     &                           lambda_stat(:,:,:,:),aa(:,:),
     &                           g_om_nu(:,:,:,:,:,:),v2_q(:,:,:),
     &                           z_red(:,:,:,:),ps3(:,:,:,:,:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           vrt_dyn_old(:,:,:,:,:,:),
     &                           vrt_stat_old(:,:,:,:),phf(:,:),
     &                           z_full(:,:,:,:),a_full(:,:,:,:)
      spin_flips=.false.
      if(theta_bse=='_GW') spin_flips=.true.
      mx2=maxel_red**2
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
c ---------------------------------------------------------------------
      call norma_wf_red
c ------ Here for zero-order vertex in full basis ------------------      
      allocate(fif_full(maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ---------- Interaction matrices in Q-orthonormal product basis ----
      allocate(w_tau(n_pbtot_red,n_pbtot_red,ndim3_tau,nqdiv_c))
      call w_tau_prepare_red(w_tau)
      if(spin_flips) call set_wnu_vrt
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
      call timel('**** P00_1Q_DYN_W : v2_q finished **')
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
      call timel('*** P00_1Q_DYN_W : all W finished **')
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
c -----------------
      if(spin_flips) then
        allocate(g_rs_mm(nrel,nfun_red,nfun_red,0:n_tau,nqdiv_c,nspin))
        allocate(g_rs_mi(nrel,nfun_red,nrel*nr_full_red,0:n_tau,nqdiv_c,
     &                   nspin))
        allocate(g_rs_im(nrel,nrel*nr_full_red,nfun_red,0:n_tau,nqdiv_c,
     &                   nspin))
        allocate(g_rs_ii(nrel,nrel*nr_red_red,nrel*nr_full_red,0:n_tau,
     &                   nqdiv_c,nspin))
        do ispin=1,nspin
          call g_real_space(g_tau(1,1,1,1,1,ispin),
     &                      g_rs_mm(1,1,1,0,1,ispin),
     &                      g_rs_mi(1,1,1,0,1,ispin),
     &                      g_rs_im(1,1,1,0,1,ispin),
     &                      g_rs_ii(1,1,1,0,1,ispin),z_red(1,1,1,ispin),
     &                      a_red(1,1,1,ispin),phf,ispin)
        enddo
        call timel('**** G_REAL_SPACE finished *********')
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
      call timel('*** P00_1Q_DYN_W : G_0 finished ****')
      allocate(k0_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2,nqdiv_c,nspin))
      allocate(k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2,nqdiv_c,nspin))
      allocate(lambda_stat(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      allocate(vrt_stat_old(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      if(spin_flips) then
        allocate(eqr(nqdiv_c))
      endif
c --------------------------------------------------------------------
      allocate(conv_table(0:iter_ladder_xi,0:7))
      conv_table=0.d0
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
      call timel('**** P00_1Q_DYN_W : loop starts ****')
      iq=q_suscept
      iqf=k_a_from_c(iq)
      if(spin_flips) call eqr_vertex_line(eqr,iq)
      ka=k_a_from_c(iq)
      np=n_pbmt_red+nplwgw_red(ka)
      allocate(pol(n_pbtot_red,n_pbtot_red,ndim3_nu))
      pol=(0.d0,0.d0)
      ps3=(0.d0,0.d0)
      do ispin=1,nspin
        do k=1,nqdiv_c
          v=pnt_c(:,k)-pnt_c(:,iq)
          call zone1_number(v,rb0_c,ndiv_c,kq)
          kq=index_k1_c(kq)    !! for K-Q
          gtild(:,1)=0.d0
          gtild(:,2)=v-pnt_c(:,kq)
          call ppm_factors_full(k,kq,iq,fif_full(1,1,1,ispin),
     &                          z_full(1,1,1,ispin),
     &                          a_full(1,1,1,ispin),
     &                          ps3(1,1,1,k,ispin),gtild,
     &                          n_low_bnd(1,ispin))
        enddo   !! over k
      enddo   !! over ispin
      lam_max=maxval(abs(ps3))
      conv_table(0,:)=lam_max/nproc_t
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        nit=iter_ladder_xi
        nom=num_omega1_adapt(i_nu)
        nomnu=num_nu2_adapt(i_nu)
        allocate(lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,nqdiv_c,nspin))
        allocate(vrt_dyn_old(nbndf_bnd,nbndf_bnd,nom,2,nqdiv_c,
     &                       nspin))
        allocate(g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,npnt_c,nspin))
        g_om_nu=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
        call g_interp_bnd(gexa,gspl_geom,gspl_asy,gasy,g_om_nu,i_nu,
     &                    nom)
c --------------------------------------------------------------------
        do ind_s2=1,ndim4_pbr(me_b+1,iq)
          s2=n4_mpi_pbr(me_b+1,iq)+ind_s2
          if(i_nu==0.and.s2==1)
     &       call timel('**** K0_BND starts *****************')
          lambda_dyn=(0.d0,0.d0)
c -------- Get zero order generalized susceptibility -------------------
          call k0_bnd(ind_nu,iq,k0_pw,g_om_nu,ps3,s2,g_x,ex,nom)
          if(i_nu==0.and.s2==1)
     &       call timel('**** K0_BND finished ***************')
c ----------------------------------------------------------------------
          k_pw=(0.d0,0.d0)
          vrt_dyn_old=(0.d0,0.d0)
          vrt_stat_old=(0.d0,0.d0)
          do it_vrt=1,nit
            if(i_nu==0.and.s2==1.and.it_vrt==1)
     &          call timel('****** Iter-n starts ***************')
            k_pw=k0_pw+k_pw
            call vertex_bnd_loop(spin_flips,iq,ind_nu,it_vrt,2,s2,nom,
     &                           nomnu,k_pw,fif0,z_red,a_red,phf,ax,
     &                           v4_r,vmi,vim,vii,w4_r,wmi,wim,wii,
     &                           v2_q,lambda_stat,lambda_dyn,nspin,0,
     &                           nrr3,indrr3)
            if(p_reducible) then
              call vertex_coulomb(iq,np,k_pw,ps3,v2_q(1,1,iq),
     &                            lambda_stat)
            endif
c --------------------------------------------------------------------
            conv_1=0.d0
            do ispin=1,nspin
              do k=1,nqdiv_c
                do i_omega=1,nom
                  do it=1,2
                    do j=1,nbndf_bnd
                      do i=1,nbndf_bnd
                        cc=lambda_dyn(i,j,i_omega,it,k,ispin)
     &                    +lambda_stat(i,j,k,ispin)
     &                    -vrt_dyn_old(i,j,i_omega,it,k,ispin)
     &                    -vrt_stat_old(i,j,k,ispin)
                        conv_1=max(conv_1,abs(cc))
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
            if(i_nu<=7) then
                conv_table(it_vrt,i_nu)=max(conv_1,
     &                                      conv_table(it_vrt,i_nu))
            endif
            vrt_dyn_old=lambda_dyn
            vrt_stat_old=lambda_stat
            do ispin=1,nspin
              if(it_vrt<nit) then
                k_pw(:,:,:,:,:,ispin)=(0.d0,0.d0)
                do k=1,nqdiv_c
c -------- Get correction to the generalized susceptibility --------
                  call k_bnd(ispin,nom,ind_nu,k,iq,
     &                       lambda_dyn(1,1,1,1,k,ispin),
     &                       k_pw(1,1,0,1,k,ispin),
     &                       lambda_stat(1,1,k,ispin),
     &                       g_om_nu(1,1,1,1,1,ispin),
     &                       g_x(1,0,1,ispin),ex(1,1,ispin))
                enddo  !! over k
              else if(it_vrt==nit) then
c -------------- Correction to the Polarizability ----------------------
                do k=1,nqdiv_c
                  call dp_bnd(nom,ind_nu,k,iq,
     &                        lambda_dyn(1,1,1,1,k,ispin),
     &                        lambda_stat(1,1,k,ispin),
     &                        g_om_nu(1,1,1,1,1,ispin),
     &                        ps3(1,1,1,k,ispin),pol(1,s2,ind_nu),ispin,
     &                        g_x(1,0,1,ispin),ex(1,1,ispin),
     &                        trn(1,i_nu),0,s2)
                enddo  !! over k
              endif
            enddo  !! over ispin
          enddo  !! over it_vrt
        enddo   !! over ind_s2
        if(i_nu==0) call timel('**** S2-loop finished **************')
        deallocate(lambda_dyn,vrt_dyn_old,g_om_nu)
      enddo   !! over ind_nu
      call timel('**** Nu-loop finished **************')
      deallocate(v4_r,w4_r,vmi,wmi,k0_pw,k_pw,trn,lambda_stat,gexa,
     &           gspl_geom,gspl_asy,gasy,g_x,ex,vim,wim,vii,wii,
     &           vrt_stat_old,fif_full,z_full,a_full,indrr3)
      if(spin_flips) deallocate(g_rs_mm,g_rs_mi,g_rs_im,g_rs_ii,
     &                          wnu_vrt,eqr)
c ---------------------------------------------------------------
      if(nproc_b/=1) call dgop(pol,2*n_pbtot_red**2*ndim3_nu,'  +',
     &                           comm_b)
      if(nspin==1.and.irel/=2) pol=2.d0*pol
      p_herm=0.d0 
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        ph=check_hermit(pol(1,1,ind_nu),np,n_pbtot_red)
        p_herm=max(p_herm,ph)
      enddo
      if(nproc_t/=1) call dgop(p_herm,1,'max',comm_t)
      if(maswrk) write(iun,*)' P00_1Q_DYN_W : P_HERM = ',p_herm
      deallocate(ps3,z_red,a_red)
      pw=(0.d0,0.d0)
      ew=(0.d0,0.d0)
      e1w=(0.d0,0.d0)
      rw=(0.d0,0.d0)
      allocate(trn(np,np))
      allocate(aa(np,np))
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(p_reducible) then
          call pw_from_pb_red(iq,np,n_pbtot_red,pol(1,1,ind_nu),
     &                        rw(i_nu),3)
	      call zgemm('n','n',np,np,np,(1.d0,0.d0),v2_q(1,1,iq),
     &               n_pbtot_red,pol(1,1,ind_nu),n_pbtot_red,
     &               (0.d0,0.d0),trn,np)
          do i=1,np
            trn(i,i)=trn(i,i)+(1.d0,0.d0)
          enddo
          call pw_from_pb_red(iq,np,np,trn,e1w(i_nu),2)
          ew(i_nu)=(1.d0,0.d0)/e1w(i_nu)
	      call zgemm('n','n',np,np,np,(1.d0,0.d0),pol(1,1,ind_nu),
     &               n_pbtot_red,trn,np,(0.d0,0.d0),aa,np)
          call pw_from_pb_red(iq,np,np,aa,pw(i_nu),3)
        else
          call pw_from_pb_red(iq,np,n_pbtot_red,pol(1,1,ind_nu),
     &                        pw(i_nu),3)
	      call zgemm('n','n',np,np,np,(1.d0,0.d0),v2_q(1,1,iq),
     &               n_pbtot_red,pol(1,1,ind_nu),n_pbtot_red,
     &               (0.d0,0.d0),trn,np)
          trn=-trn
          do i=1,np
            trn(i,i)=trn(i,i)+(1.d0,0.d0)
          enddo
          call invers_z(np,trn,np)
          call pw_from_pb_red(iq,np,np,trn,e1w(i_nu),2)
          ew(i_nu)=(1.d0,0.d0)/e1w(i_nu)
	      call zgemm('n','n',np,np,np,(1.d0,0.d0),pol(1,1,ind_nu),
     &               n_pbtot_red,trn,np,(0.d0,0.d0),aa,np)
          call pw_from_pb_red(iq,np,np,aa,rw(i_nu),3)
        endif
      enddo
      deallocate(v2_q,trn,aa)
      if(nproc_t/=1) then
        call dgop(pw,2*(n_nu+1),'  +',comm_t)
        call dgop(rw,2*(n_nu+1),'  +',comm_t)
        call dgop(ew,2*(n_nu+1),'  +',comm_t)
        call dgop(e1w,2*(n_nu+1),'  +',comm_t)
      endif
c ------------------------------------------------------------------
      idum=(iter_ladder_xi+1)*8
      if(nproc_b/=1) call dgop(conv_table,idum,'max',comm_b)
      if(nproc_t/=1) call dgop(conv_table,idum,'  +',comm_t)
      call output_bs_conv_1q(conv_table,iter_ladder_xi)
      deallocate(pol,conv_table,phf,nind_fif,ind_fif,fif0,ax)
      call timel('**** P00_1Q_DYN_W finished *********')
      end
