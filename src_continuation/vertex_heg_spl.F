      subroutine vertex_heg_spl(p_need,sig_need,iter_lad,
     &                          spin_flips)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical, intent(in) :: p_need,sig_need,spin_flips
      integer, intent(in) :: iter_lad
      logical :: neg
      integer :: ind_nu,i_nu,i,nom,i_omega,nom_nu,i_tau,ind_tau,it,kma,
     &           ind_q,iq,k,ig,k0,i_wi,iter_vrt,ispin,iqc,k1
      real*8 :: w1,w3,conv_1
      complex*16 :: cc
      real*8, allocatable :: wspl_geom(:,:,:),wspl_asy(:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:),gasy(:,:,:),
     &                       gspl_asy(:,:,:,:,:),pol(:,:),gexa(:,:,:,:),
     &                       s(:,:,:),sa1(:,:,:),conv_table(:,:,:),
     &                       conv_cur_table(:,:,:),s_dyn(:,:,:),
     &                       w_q_nuc_heg(:,:),
     &                       wexa(:,:),wasy(:)
      complex*16, allocatable :: g_om(:,:,:,:),g_om_x(:,:,:,:),
     &                           lambda_stat_nu(:,:,:,:),
     &                           lambda_stat_nu1(:,:,:,:),
     &                           lambda_stat_nu_sig2(:,:,:,:),
     &                           k_pw(:,:,:),vrt_out(:,:),
     &                           aw(:,:),bw(:,:),
     &                           lambda_dyn(:,:,:),
     &                           lambda_dyn_2(:,:,:),
     &                           lambda_stat_tau(:,:,:,:,:),
     &                           vrt_k_q(:,:),vrt_nu_q(:,:),
     &                           vrt_om_nu(:,:),vrt_om_q(:,:),
     &                           vrt_om_k(:,:),vrt_nu_k(:,:)
      kma=min(3,nc_line)
      if(sig_need) then
        allocate(s(nq_sph_c,2,ndim3_tau))
        allocate(s_dyn(nq_sph_c,0:n_tau,nspin))
        allocate(sa1(nq_sph_c,2,ndim3_tau))
      endif
      allocate(lambda_stat_nu_sig2(nrr_div_c,nspin,ndim3_nu,
     &                             ndim_k_red(me_k+1)))
c ----------------------------------------------------------------
      allocate(conv_table(2:iter_lad,0:7,kma))
      allocate(conv_cur_table(2:iter_ladder_p,0:7,kma))
      allocate(lambda_stat_nu(nrr_div_c,nspin,ndim3_nu,
     &	                    ndim_k_red(me_k+1)))
      allocate(vrt_out(nc_line,0:n_nu))
      allocate(vrt_k_q(nc_line,nc_line))
      allocate(vrt_nu_q(0:n_nu_exa,nc_line))
      allocate(vrt_om_nu(0:n_omega_exa,0:n_nu_exa))
      allocate(vrt_om_q(0:n_omega_exa,nc_line))
      allocate(vrt_om_k(0:n_omega_exa,nc_line))
      allocate(vrt_nu_k(0:n_nu_exa,nc_line))
c --------------------------------------------------------------------
      allocate(k_pw(nrr_div_c,0:n_tau,nspin))
      allocate(lambda_stat_tau(nrr_div_c,2,ndim3_tau,
     &                         ndim_k_red(me_k+1),nspin))
c ------------------------------------------------------------------
      allocate(gexa(0:n_omega_exa,2,nstar_c,nspin))
      allocate(gspl_geom(0:n_omega_geom+2,4,2,nstar_c,nspin))
      allocate(gspl_asy(n_omega_asy+1,4,2,nstar_c,nspin))
      allocate(gasy(2,nstar_c,nspin))
      allocate(g0(0:n_omega,2))
      do ispin=1,nspin
        do i=1,nstar_c
          do i_omega=0,n_omega
            cc=g_k_omega_heg_c(i_omega,i,ispin)
     &	    -gx_k_omega_heg_c(i_omega,i,ispin)
            g0(i_omega,1)=real(cc)
            g0(i_omega,2)=imag(cc)
          enddo
          do it=1,2
            call spline_inhmg(w_omega(n_omega_exa-1),
     &                        g0(n_omega_exa-1,it),
     &                        gspl_geom(0,1,it,i,ispin),
     &                        gspl_geom(0,2,it,i,ispin),
     &                        gspl_geom(0,3,it,i,ispin),
     &                        gspl_geom(0,4,it,i,ispin),n_omega_geom+3,
     &                        0,0.d0,0.d0)
            call spline_inhmg(xm_omega(0),
     &                        g0(n_omega_exa+n_omega_geom,it),
     &                        gspl_asy(1,1,it,i,ispin),
     &                        gspl_asy(1,2,it,i,ispin),
     &                        gspl_asy(1,3,it,i,ispin),
     &                        gspl_asy(1,4,it,i,ispin),n_omega_asy+1,
     &                        0,0.d0,0.d0)
            gasy(it,i,ispin)=g0(n_omega,it)
            gexa(:,it,i,ispin)=g0(0:n_omega_exa,it)
          enddo
        enddo
      enddo
      deallocate(g0)
c -----------------------------------------------------------------
      allocate(wexa(0:n_nu_exa,nstar_c))
      allocate(wspl_geom(0:n_nu_geom+2,4,nstar_c))
      allocate(wspl_asy(n_nu_asy+1,4,nstar_c))
      allocate(wasy(nstar_c))
      if(spin_flips.or.sig_need) then
        do i=1,nstar_c
          call spline_inhmg(w_nu(n_nu_exa-1),
     &                      w_q_nu_heg_c(n_nu_exa-1,i),
     &                      wspl_geom(0,1,i),wspl_geom(0,2,i),
     &                      wspl_geom(0,3,i),wspl_geom(0,4,i),
     &                      n_nu_geom+3,0,0.d0,0.d0)
          call spline_inhmg(xm_nu(0),
     &                      w_q_nu_heg_c(n_nu_exa+n_nu_geom,i),
     &                      wspl_asy(1,1,i),wspl_asy(1,2,i),
     &                      wspl_asy(1,3,i),wspl_asy(1,4,i),n_nu_asy+1,
     &                      0,0.d0,0.d0)
          wexa(:,i)=w_q_nu_heg_c(0:n_nu_exa,i)
          wasy(i)=w_q_nu_heg_c(n_nu,i)
        enddo
      endif
c --------------------------------------------------------------------
      if(p_need) then
        allocate(pol(nc_line,0:n_nu))
        pol=0.d0
        vrt_out=(0.d0,0.d0)
        vrt_k_q=(0.d0,0.d0)
        vrt_nu_q=(0.d0,0.d0)
        vrt_om_nu=(0.d0,0.d0)
        vrt_om_q=(0.d0,0.d0)
        vrt_om_k=(0.d0,0.d0)
        vrt_nu_k=(0.d0,0.d0)
      endif
c --------------------------------------------------------------------
      if(sig_need) then
        s_dyn=0.d0
      endif
      if(iter_lad>1) conv_table=0.d0
      if(iter_ladder_p>1) conv_cur_table=0.d0
      conv_ladder=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        nom=num_omega1_adapt(i_nu)
c --------------------------------------------------------------------
        w1=w_nu(i_nu)-pi/betta_t
        i_wi=1
        neg=.false.
        if(w1>0.d0) then
          neg=.true.
          do i_omega=1,nom
            if(abs(w_om_adapt_nu(i_omega,i_nu)-w1)<1.d-7) i_wi=i_omega
          enddo
        endif
c --------------------------------------------------------------------
        nom_nu=num_nu2_adapt(i_nu)
        allocate(lambda_dyn(nrr_div_c,nom,nspin))
        allocate(lambda_dyn_2(nrr_div_c,nom,nspin))
        allocate(g_om(nstar_c,nom,2,nspin))
        allocate(g_om_x(nstar_c,nom,2,nspin))
c ------ G - interpolation -----------------------------------
        do ispin=1,nspin
          do i_omega=1,nom
            w1=w_om_adapt_nu(i_omega,i_nu)
            call interp_g_spl(gexa(0,1,1,ispin),
     &                        gspl_geom(0,1,1,1,ispin),
     &                        gspl_asy(1,1,1,1,ispin),gasy(1,1,ispin),
     &                        g_om(1,i_omega,1,ispin),w1,nstar_c,3)
            do i=1,nstar_c
              g_om_x(i,i_omega,1,ispin)=(1.d0,0.d0)
     &                       /dcmplx(chem_pot-e_star_xc(i,ispin),w1)
            enddo
            w1=w_om_adapt_nu(i_omega,i_nu)-w_nu(i_nu)
            call interp_g_spl(gexa(0,1,1,ispin),
     &                        gspl_geom(0,1,1,1,ispin),
     &                        gspl_asy(1,1,1,1,ispin),gasy(1,1,ispin),
     &                        g_om(1,i_omega,2,ispin),w1,nstar_c,3)
            do i=1,nstar_c
              g_om_x(i,i_omega,2,ispin)=(1.d0,0.d0)
     &                        /dcmplx(chem_pot-e_star_xc(i,ispin),w1)
            enddo
          enddo
        enddo
        g_om=g_om+g_om_x
c ----------------------------------------------------------------
        do ind_q=1,ndim_k_red(me_k+1)
          iq=n_mpi_k_red(me_k+1)+ind_q
          iqc=k_npnt_in_line_c(iq)
          iter_vrt=iter_lad
! ---------------------------------------------------------------- 
          call vertex_heg_a_spl(iter_vrt,iq,ind_nu,nom,lambda_dyn,
     &	                      lambda_dyn_2,
     &                          lambda_stat_nu(1,1,ind_nu,ind_q),
     &                          lambda_stat_nu_sig2(1,1,ind_nu,ind_q),
     &                          g_om,g_om_x,wexa,wspl_geom,wspl_asy,
     &                          wasy,conv_table,conv_1,spin_flips,
     &                          iter_lad,q_irr_c(1,iq),kma)
          if(p_need.and.iqc/=0) conv_ladder(iqc,i_nu)=conv_1
          if(sig_need) then
            call interp_w_nu_spl(wexa(0,iq),wspl_geom(0,1,iq),
     &                           wspl_asy(1,1,iq),wasy(iq),
     &                           w3,w_nu(i_nu),1)
            allocate(aw(nq_sph_c,0:n_tau))
            allocate(bw(nq_sph_c,0:n_tau))
            do ispin=1,nspin
              aw=(0.d0,0.d0)
              bw=(0.d0,0.d0)
              call sigma_g_lambda_heg_0(iq,nom,ind_nu,
     &                                  lambda_dyn(1,1,ispin),
     &                                  g_om(1,1,1,ispin),w3,aw,bw)
              if(iter_sigma_gwg>1) then
                call sigma_g_lambda_heg_0(iq,nom,ind_nu,
     &                                    lambda_dyn_2(1,1,ispin),
     &                                    g_om(1,1,1,ispin),
     &                                    v_q_heg_c(iq),aw,bw)
              endif
              call sigma_vrt_dyn_heg(aw,bw,s_dyn(1,0,ispin),ind_nu)
            enddo
            deallocate(aw,bw)
          endif
          if(p_need.and.iqc/=0) then
c -------------- Correction to the Polarizability ----------------------
            do ispin=1,nspin
              call k_heg(ispin,q_irr_c(1,iq),nom,ind_nu,
     &                   lambda_dyn(1,1,ispin),k_pw(1,0,ispin),
     &                   lambda_stat_nu(1,ispin,ind_nu,ind_q),
     &                   g_om(1,1,1,ispin),g_om_x(1,1,1,ispin))
              pol(iqc,i_nu)=pol(iqc,i_nu)-2.d0*real(k_pw(1,0,ispin))
            enddo
c ----- Prepare Vertex for output ----------------------------------
            call heg_vrt_for_out(vrt_out,vrt_k_q,vrt_om_nu,vrt_om_q,
     &                           vrt_om_k,vrt_nu_k,vrt_nu_q,ind_q,
     &                           ind_nu,nom,neg,i_wi,
     &                           lambda_stat_nu(1,1,ind_nu,ind_q),
     &                           lambda_dyn)
          endif
        enddo   !! over ind_q
        deallocate(lambda_dyn,g_om,g_om_x,lambda_dyn_2)
      enddo   !! over ind_nu
      if(p_need.and.nspin==1) pol=2.d0*pol
      if(nproc_k/=1) then
        if(sig_need) then
          call dgop(s_dyn,nq_sph_c*(n_tau+1)*nspin,'  +',comm_k)
        endif
        if(p_need) then
          call dgop(pol,nc_line*(1+n_nu),'  +',comm_k)
          call dgop(conv_ladder,nc_line*(1+n_nu),'  +',comm_k)
          call dgop(vrt_out,2*nc_line*(1+n_nu),'  +',comm_k)
          call dgop(vrt_k_q,2*nc_line**2,'  +',comm_k)
          call dgop(vrt_nu_q,2*nc_line*(1+n_nu_exa),'  +',comm_k)
          call dgop(vrt_om_nu,2*(1+n_omega_exa)*(1+n_nu_exa),'  +',
     &              comm_k)
          call dgop(vrt_om_q,2*nc_line*(1+n_omega_exa),'  +',comm_k)
          call dgop(vrt_om_k,2*nc_line*(1+n_omega_exa),'  +',comm_k)
          call dgop(vrt_nu_k,2*nc_line*(1+n_nu_exa),'  +',comm_k)
          if(iter_lad>1) call dgop(conv_table,(iter_lad-1)*8*kma,'  +',
     &	                         comm_k)
        endif
      endif
      if(nproc_t/=1) then
        if(p_need) then
          call dgop(pol,nc_line*(1+n_nu),'  +',comm_t)
          call dgop(conv_ladder,nc_line*(1+n_nu),'  +',comm_t)
          call dgop(vrt_out,2*nc_line*(1+n_nu),'  +',comm_t)
          call dgop(vrt_k_q,2*nc_line**2,'  +',comm_t)
          call dgop(vrt_nu_q,2*nc_line*(1+n_nu_exa),'  +',comm_t)
          call dgop(vrt_om_nu,2*(1+n_omega_exa)*(1+n_nu_exa),'  +',
     &              comm_t)
          call dgop(vrt_om_q,2*nc_line*(1+n_omega_exa),'  +',comm_t)
          call dgop(vrt_om_k,2*nc_line*(1+n_omega_exa),'  +',comm_t)
          call dgop(vrt_nu_k,2*nc_line*(1+n_nu_exa),'  +',comm_t)
          if(iter_lad>1) call dgop(conv_table,(iter_ladder_p-1)*8*kma,
     &	                         '  +',comm_t)
        endif
        if(sig_need) then
        call dgop(s_dyn,nq_sph_c*(n_tau+1)*nspin,'  +',comm_t)
        endif
      endif
      if(maswrk.and.p_need.and.iter_lad>1) then
        do k=1,kma
          write(iun,*)' Convergency of the ladder sequence, Q =',k
          write(iun,'(a5,8(i3,6x))')'It\Nu',0,1,2,3,4,5,6,7
          do it=2,iter_lad
            write(iun,'(i3,8e9.2)')it,(conv_table(it,i_nu,k),i_nu=0,7)
          enddo
        enddo
      endif
c ------------------------------------------------------------------
      if(p_need) then
        call output_heg_vrt_spl(vrt_out,vrt_k_q,vrt_nu_q,vrt_om_nu,
     &                          vrt_om_q,vrt_om_k,vrt_nu_k)
c ------- Interpolation of P onto fine mesh -----------------------
        do i_nu=0,n_nu
          do i=1,nc_line
            p_q_nu_vertex_heg(i_nu,i)=pol(i,i_nu)
          enddo
        enddo
      endif
      allocate(lambda_stat_nu1(nrr_div_c,ndim3_nu,ndim_k_red(me_k+1),
     &                         nspin))
      call vertex_heg_stat_spl(lambda_stat_tau,
     &	                     lambda_stat_nu_sig2,
     &	                     lambda_stat_nu1)
      if(sig_need) then
	  allocate(w_q_nuc_heg(ndim3_nu,nstar_c))
	  do i=1,nstar_c
	    do ind_nu=1,ndim3_nu
	      i_nu=me_t*ndim3_nu+ind_nu-1
            w1=w_nu(i_nu)
            call interp_w_nu_spl(wexa(0,i),wspl_geom(0,1,i),
     &                           wspl_asy(1,1,i),wasy(i),
     &                           w_q_nuc_heg(ind_nu,i),w1,1)
          enddo
        enddo
        do ispin=1,nspin
c ----- Semidynamic Sigma A of the first only order -----------------
          call sig_heg_semi(ispin,lambda_stat_nu1(1,1,1,ispin),sa1)
c ------------- The rest of Sigma ---------------------------------          
          call sigma_g_lambda_heg_spl(ispin,
     &	                            lambda_stat_tau(1,1,1,1,ispin),
     &                                lambda_stat_nu_sig2,s,w_q_nuc_heg)
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              if(it==2.and.i_tau==n_tau/2) cycle
              s(:,it,ind_tau)=s(:,it,ind_tau)+sa1(:,it,ind_tau)
     &                       +s_dyn(:,i_tau,ispin)
c ------------ Symmetrization -------------------------------------
              do k=1,nc_line
                k1=kline_in_npnt_c(k)
                do ig=1,ngroup
                  k0=ip_k_c(k1,ig)
                  sig_k_tau_vertex_heg(i_tau,k,ispin)=
     &              sig_k_tau_vertex_heg(i_tau,k,ispin)+s(k0,it,ind_tau)
                enddo
              enddo
            enddo
          enddo
        enddo
        deallocate(w_q_nuc_heg)
        if(nproc_t/=1) call dgop(sig_k_tau_vertex_heg,
     & 	                         nc_line*(1+n_tau)*nspin,'  +',
     &                             comm_t)
        sig_k_tau_vertex_heg=sig_k_tau_vertex_heg/ngroup
      endif
      deallocate(k_pw,gspl_geom,gspl_asy,lambda_stat_nu,gasy,
     &           lambda_stat_tau,gexa,wexa,wasy,
     &           lambda_stat_nu1,conv_table,wspl_geom,wspl_asy,
     &           lambda_stat_nu_sig2,vrt_out,vrt_k_q,vrt_nu_q,
     &           vrt_om_nu,vrt_om_q,vrt_om_k,vrt_nu_k)
      if(p_need) deallocate(pol)
      if(sig_need) deallocate(s,sa1,s_dyn)
      if(p_need) then
        call output_ladder_conv(0)
      endif
      call timel('**** VERTEX_HEG_SPL finished *******')
      end