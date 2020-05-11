      subroutine vertex_heg_p_only_spl(iter_lad,spin_flips)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical, intent(in) :: spin_flips
      integer, intent(in) :: iter_lad
      logical :: neg
      integer :: ind_nu,i_nu,i,nom,i_omega,nom_nu,it,kma,
     &           ind_q,iq,k,i_wi,iter_vrt,ispin,iqc
      real*8 :: w1,conv_1
      complex*16 :: cc
      integer, allocatable :: ndimc_k_line(:),nc_mpi_k_line(:)
      real*8, allocatable :: wspl_geom(:,:,:),wspl_asy(:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:),gasy(:,:,:),
     &                       gspl_asy(:,:,:,:,:),pol(:,:),gexa(:,:,:,:),
     &                       conv_table(:,:,:),
     &                       conv_cur_table(:,:,:),
     &                       wexa(:,:),wasy(:)
      complex*16, allocatable :: g_om(:,:,:,:),g_om_x(:,:,:,:),
     &                           lambda_stat_nu(:,:),
     &                           lambda_stat_nu_sig2(:,:),
     &                           k_pw(:,:,:),vrt_out(:,:),
     &                           lambda_dyn(:,:,:),
     &                           lambda_dyn_2(:,:,:),
     &                           vrt_k_q(:,:),vrt_nu_q(:,:),
     &                           vrt_om_nu(:,:),vrt_om_q(:,:),
     &                           vrt_om_k(:,:),vrt_nu_k(:,:)
      kma=min(3,nc_line)
      allocate(ndimc_k_line(nproc_k))
      allocate(nc_mpi_k_line(nproc_k))
      call size_shift_par(nc_line,nproc_k,ndimc_k_line,nc_mpi_k_line)
      allocate(lambda_stat_nu_sig2(nrr_div_c,nspin))
c ----------------------------------------------------------------
      allocate(conv_table(2:iter_lad,0:7,kma))
      allocate(conv_cur_table(2:iter_ladder_p,0:7,kma))
      allocate(lambda_stat_nu(nrr_div_c,nspin))
      allocate(vrt_out(nc_line,0:n_nu))
      allocate(vrt_k_q(nc_line,nc_line))
      allocate(vrt_nu_q(0:n_nu_exa,nc_line))
      allocate(vrt_om_nu(0:n_omega_exa,0:n_nu_exa))
      allocate(vrt_om_q(0:n_omega_exa,nc_line))
      allocate(vrt_om_k(0:n_omega_exa,nc_line))
      allocate(vrt_nu_k(0:n_nu_exa,nc_line))
c --------------------------------------------------------------------
      allocate(k_pw(nrr_div_c,0:n_tau,nspin))
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
      if(spin_flips) then
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
      allocate(pol(nc_line,0:n_nu))
      pol=0.d0
      vrt_out=(0.d0,0.d0)
      vrt_k_q=(0.d0,0.d0)
      vrt_nu_q=(0.d0,0.d0)
      vrt_om_nu=(0.d0,0.d0)
      vrt_om_q=(0.d0,0.d0)
      vrt_om_k=(0.d0,0.d0)
      vrt_nu_k=(0.d0,0.d0)
c --------------------------------------------------------------------
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
        do ind_q=1,ndimc_k_line(me_k+1)
          iqc=nc_mpi_k_line(me_k+1)+ind_q
          iq=kline_in_npnt_c(iqc)
          iter_vrt=iter_lad
! ---------------------------------------------------------------- 
          call vertex_heg_a_spl(iter_vrt,iq,ind_nu,nom,lambda_dyn,
     &	                      lambda_dyn_2,lambda_stat_nu,
     &                          lambda_stat_nu_sig2,
     &                          g_om,g_om_x,wexa,wspl_geom,wspl_asy,
     &                          wasy,conv_table,conv_1,spin_flips,
     &                          iter_lad,q_irr_c(1,iq),kma)
          conv_ladder(iqc,i_nu)=conv_1
c -------------- Correction to the Polarizability ----------------------
          do ispin=1,nspin
            call k_heg(ispin,q_irr_c(1,iq),nom,ind_nu,
     &                 lambda_dyn(1,1,ispin),k_pw(1,0,ispin),
     &                 lambda_stat_nu(1,ispin),
     &                 g_om(1,1,1,ispin),g_om_x(1,1,1,ispin))
            pol(iqc,i_nu)=pol(iqc,i_nu)-2.d0*real(k_pw(1,0,ispin))
          enddo
c ----- Prepare Vertex for output ----------------------------------
          call heg_vrt_for_out_line(vrt_out,vrt_k_q,vrt_om_nu,vrt_om_q,
     &                              vrt_om_k,vrt_nu_k,vrt_nu_q,iqc,
     &                              ind_nu,nom,neg,i_wi,lambda_stat_nu,
     &                              lambda_dyn)
        enddo   !! over ind_q
        deallocate(lambda_dyn,g_om,g_om_x,lambda_dyn_2)
      enddo   !! over ind_nu
      if(nspin==1) pol=2.d0*pol
      if(nproc_k/=1) then
        call dgop(pol,nc_line*(1+n_nu),'  +',comm_k)
        call dgop(conv_ladder,nc_line*(1+n_nu),'  +',comm_k)
        call dgop(vrt_out,2*nc_line*(1+n_nu),'  +',comm_k)
        call dgop(vrt_k_q,2*nc_line**2,'  +',comm_k)
        call dgop(vrt_nu_q,2*nc_line*(1+n_nu_exa),'  +',comm_k)
        call dgop(vrt_om_nu,2*(1+n_omega_exa)*(1+n_nu_exa),'  +',
     &            comm_k)
        call dgop(vrt_om_q,2*nc_line*(1+n_omega_exa),'  +',comm_k)
        call dgop(vrt_om_k,2*nc_line*(1+n_omega_exa),'  +',comm_k)
        call dgop(vrt_nu_k,2*nc_line*(1+n_nu_exa),'  +',comm_k)
        if(iter_lad>1) call dgop(conv_table,(iter_lad-1)*8*kma,'  +',
     &	                       comm_k)
      endif
      if(nproc_t/=1) then
        call dgop(pol,nc_line*(1+n_nu),'  +',comm_t)
        call dgop(conv_ladder,nc_line*(1+n_nu),'  +',comm_t)
        call dgop(vrt_out,2*nc_line*(1+n_nu),'  +',comm_t)
        call dgop(vrt_k_q,2*nc_line**2,'  +',comm_t)
        call dgop(vrt_nu_q,2*nc_line*(1+n_nu_exa),'  +',comm_t)
        call dgop(vrt_om_nu,2*(1+n_omega_exa)*(1+n_nu_exa),'  +',
     &            comm_t)
        call dgop(vrt_om_q,2*nc_line*(1+n_omega_exa),'  +',comm_t)
        call dgop(vrt_om_k,2*nc_line*(1+n_omega_exa),'  +',comm_t)
        call dgop(vrt_nu_k,2*nc_line*(1+n_nu_exa),'  +',comm_t)
        if(iter_lad>1) call dgop(conv_table,(iter_lad-1)*8*kma,
     &	                       '  +',comm_t)
      endif
      if(maswrk.and.iter_lad>1) then
        do k=1,kma
          write(iun,*)' Convergency of the ladder sequence, Q =',k
          write(iun,'(a5,8(i3,6x))')'It\Nu',0,1,2,3,4,5,6,7
          do it=2,iter_lad
            write(iun,'(i3,8e9.2)')it,(conv_table(it,i_nu,k),i_nu=0,7)
          enddo
        enddo
      endif
c ------------------------------------------------------------------
      call output_heg_vrt_spl(vrt_out,vrt_k_q,vrt_nu_q,vrt_om_nu,
     &                        vrt_om_q,vrt_om_k,vrt_nu_k)
c ------- Interpolation of P onto fine mesh -----------------------
      do i_nu=0,n_nu
        do i=1,nc_line
          p_q_nu_vertex_heg(i_nu,i)=pol(i,i_nu)
        enddo
      enddo
      deallocate(k_pw,gspl_geom,gspl_asy,lambda_stat_nu,gasy,
     &           gexa,wexa,wasy,ndimc_k_line,nc_mpi_k_line,
     &           conv_table,wspl_geom,wspl_asy,
     &           lambda_stat_nu_sig2,vrt_out,vrt_k_q,vrt_nu_q,
     &           vrt_om_nu,vrt_om_q,vrt_om_k,vrt_nu_k)
      deallocate(pol)
      call output_ladder_conv(0)
      call timel('**** VERTEX_HEG_P_ONLY finished ****')
      end