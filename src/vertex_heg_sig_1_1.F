      subroutine vertex_heg_sig_1_1
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i_nu,i,nom,i_omega,i_tau,ind_tau,it,ispin,ind_omega,
     &           i1_omega,nomnu,ind_k,numtn,k1l,k1
      real*8 :: w1,pib
      complex*16 :: cc,g_x
      integer, allocatable :: ndimc_k_line(:),nc_mpi_k_line(:)
      real*8, allocatable :: wspl_geom(:,:,:),wspl_asy(:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:),gasy(:,:,:),
     &                       gspl_asy(:,:,:,:,:),gexa(:,:,:,:),
     &                       wexa(:,:),wasy(:),sig_stat(:,:,:,:),
     &                       sig_semi(:,:,:,:),
     &                       sig_dyn(:,:,:,:),sigt(:,:)
      complex*16, allocatable :: k_pw(:,:,:),t_pw(:,:,:),
     &                           g_om_om(:,:,:),g_nu_om(:,:,:),
     &                           sig_dyn_om(:,:,:)
      pib=pi/betta_t
      allocate(ndimc_k_line(nproc_k))
      allocate(nc_mpi_k_line(nproc_k))
      call size_shift_par(nc_line,nproc_k,ndimc_k_line,nc_mpi_k_line)
      allocate(k_pw(nrr_div_c,0:n_tau,2))
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
      do i=1,nstar_c
        call spline_inhmg(w_nu(n_nu_exa-1),
     &                    w_q_nu_heg_c(n_nu_exa-1,i),
     &                    wspl_geom(0,1,i),wspl_geom(0,2,i),
     &                    wspl_geom(0,3,i),wspl_geom(0,4,i),
     &                    n_nu_geom+3,0,0.d0,0.d0)
        call spline_inhmg(xm_nu(0),
     &                    w_q_nu_heg_c(n_nu_exa+n_nu_geom,i),
     &                    wspl_asy(1,1,i),wspl_asy(1,2,i),
     &                    wspl_asy(1,3,i),wspl_asy(1,4,i),n_nu_asy+1,
     &                    0,0.d0,0.d0)
        wexa(:,i)=w_q_nu_heg_c(0:n_nu_exa,i)
        wasy(i)=w_q_nu_heg_c(n_nu,i)
      enddo
c --------------------------------------------------------------------
      allocate(sig_dyn_om(ndim3_omega,ndimc_k_line(me_k+1),nspin))
      sig_dyn_om=(0.d0,0.d0)
      allocate(sig_stat(2,ndim3_tau,ndimc_k_line(me_k+1),nspin))
      sig_stat=0.d0
      allocate(sig_semi(2,ndim3_tau,ndimc_k_line(me_k+1),nspin))
      sig_semi=0.d0
      do ispin=1,nspin
        do ind_omega=1,ndim3_omega
          i1_omega=me_t*ndim3_omega+ind_omega-1
          nom=num_omega1_omega(i1_omega)
          nomnu=num_nu_omega(i1_omega)
          numtn=max(nomnu,n_tau+1)
          allocate(t_pw(nrr_div_c,numtn,2))
          allocate(g_om_om(nstar_c,2,nom))
          g_om_om=(0.d0,0.d0)
          allocate(g_nu_om(nstar_c,2,nomnu))
          g_nu_om=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
          do i_omega=1,nom
            w1=w_om_adapt_omega(i_omega,i1_omega)
            call interp_g_spl(gexa(0,1,1,ispin),
     &                        gspl_geom(0,1,1,1,ispin),
     &                        gspl_asy(1,1,1,1,ispin),gasy(1,1,ispin),
     &                        g_om_om(1,2,i_omega),w1,nstar_c,3)
            do i=1,nstar_c
              g_x=(1.d0,0.d0)/dcmplx(chem_pot-e_star_xc(i,ispin),w1)
              g_om_om(i,2,i_omega)=g_om_om(i,2,i_omega)+g_x
            enddo
            w1=-w_om_adapt_omega(i_omega,i1_omega)+w_omega(i1_omega)-pib
            call interp_g_spl(gexa(0,1,1,ispin),
     &                        gspl_geom(0,1,1,1,ispin),
     &                        gspl_asy(1,1,1,1,ispin),gasy(1,1,ispin),
     &                        g_om_om(1,1,i_omega),w1,nstar_c,3)
            do i=1,nstar_c
              g_x=(1.d0,0.d0)/dcmplx(chem_pot-e_star_xc(i,ispin),w1)
              g_om_om(i,1,i_omega)=g_om_om(i,1,i_omega)+g_x
            enddo
          enddo
	    do i_nu=1,nomnu
            w1=w_omega(i1_omega)-w_nu_adapt_omega(i_nu,i1_omega)
            call interp_g_spl(gexa(0,1,1,ispin),
     &                        gspl_geom(0,1,1,1,ispin),
     &                        gspl_asy(1,1,1,1,ispin),gasy(1,1,ispin),
     &                        g_nu_om(1,2,i_nu),w1,nstar_c,3)
            do i=1,nstar_c
              g_x=(1.d0,0.d0)/dcmplx(chem_pot-e_star_xc(i,ispin),w1)
              g_nu_om(i,2,i_nu)=g_nu_om(i,2,i_nu)+g_x
            enddo
            w1=w_nu_adapt_omega(i_nu,i1_omega)+pib
            call interp_g_spl(gexa(0,1,1,ispin),
     &                        gspl_geom(0,1,1,1,ispin),
     &                        gspl_asy(1,1,1,1,ispin),gasy(1,1,ispin),
     &                        g_nu_om(1,1,i_nu),w1,nstar_c,3)
            do i=1,nstar_c
              g_x=(1.d0,0.d0)/dcmplx(chem_pot-e_star_xc(i,ispin),w1)
              g_nu_om(i,1,i_nu)=g_nu_om(i,1,i_nu)+g_x
            enddo
          enddo  !! over i_nu
c ----------------------------------------------------------------
          do ind_k=1,ndimc_k_line(me_k+1)
            k1l=nc_mpi_k_line(me_k+1)+ind_k
            k1=kline_in_npnt_c(k1l)
            call k_sig_heg(ind_omega,k1,k_pw,g_om_om,nom,wexa,
     &                     wspl_geom,wspl_asy,wasy)
            call tsig_heg(ispin,k_pw,t_pw,numtn,ind_omega,nomnu)
            call tsig_end_dyn_heg(nomnu,ind_omega,k1,t_pw,g_nu_om,
     &                            sig_dyn_om(ind_omega,ind_k,ispin),
     &                            wexa,wspl_geom,wspl_asy,wasy,numtn)
          enddo   !! over ind_k
          deallocate(t_pw,g_om_om,g_nu_om)
        enddo  !! over ind_omega
c ------- Static and Semidynamic parts -------------------------------
        allocate(t_pw(nrr_div_c,2,ndim3_tau))
        do ind_k=1,ndimc_k_line(me_k+1)
          k1l=nc_mpi_k_line(me_k+1)+ind_k
          k1=kline_in_npnt_c(k1l)
          do ind_tau=1,ndim3_tau
            do it=1,2
              call sig_1_stat_heg(it,ind_tau,t_pw(1,it,ind_tau),k1,
     &                            ispin,
     &                            sig_stat(it,ind_tau,ind_k,ispin))
            enddo  !! over it
          enddo   !!! over ind_tau
          call sig_1_semi_heg(t_pw,k1,ispin,sig_semi(1,1,ind_k,ispin))
        enddo   !! over ind_k
        deallocate(t_pw)
      enddo   !! over ispin
      deallocate(wexa,wspl_geom,wspl_asy,wasy,gexa,gspl_geom,gspl_asy,
     &           gasy,k_pw)
c ------------------------------------------------------------------
      sig_k_tau_vertex_heg=0.d0
      allocate(sig_dyn(2,ndim3_tau,ndimc_k_line(me_k+1),nspin))
      sig_dyn=0.d0
      allocate(sigt(2,ndim3_tau))
      do ispin=1,nspin
        do ind_k=1,ndimc_k_line(me_k+1)
          k1l=nc_mpi_k_line(me_k+1)+ind_k
          call omega_to_tau_baa(sig_dyn_om(1,ind_k,ispin),1,sigt,1,1)
          do ind_tau=1,ndim3_tau
            sig_dyn(1,ind_tau,ind_k,ispin)=0.5d0*(sigt(1,ind_tau)
     &                                           +sigt(2,ind_tau))
            sig_dyn(2,ind_tau,ind_k,ispin)=0.5d0*(sigt(2,ind_tau)
     &                                           -sigt(1,ind_tau))
            i_tau=me_t*ndim3_tau+ind_tau-1
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              sig_semi(it,ind_tau,ind_k,ispin)=
     &          2.d0*sig_semi(it,ind_tau,ind_k,ispin)
              sig_k_tau_vertex_heg(i_tau,k1l,ispin)=
     &              sig_stat(it,ind_tau,ind_k,ispin)
     &             +sig_semi(it,ind_tau,ind_k,ispin)
     &             +sig_dyn(it,ind_tau,ind_k,ispin)
            enddo
          enddo
        enddo
      enddo
      deallocate(sigt)
      if(nproc_t/=1) call dgop(sig_k_tau_vertex_heg,
     & 	                       nc_line*(1+n_tau)*nspin,'  +',
     &                           comm_t)
      if(nproc_k/=1) call dgop(sig_k_tau_vertex_heg,
     & 	                     nc_line*(1+n_tau)*nspin,'  +',
     &                         comm_k)
      deallocate(sig_stat,sig_semi,sig_dyn,ndimc_k_line,nc_mpi_k_line)
      call timel('**** VERTEX_HEG_SIG_1_1 finished ***')
      end
