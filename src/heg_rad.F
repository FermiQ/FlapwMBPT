      subroutine heg_rad
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical :: last
      integer :: i_tau,i_omega,i,it0,ispin,iq
      real*8 :: de,green0,deltae,w,pi2a,qs,vol0,pi2,pi8,hh,q
      complex*16 :: gxx,c1,gxs
      real*8, allocatable :: pl(:),spl(:,:)
      iter_sigma_gwg=min(iter_sigma_gwg,iter_ladder_p)
      pi2=pi+pi
      pi2a=pi2/par
      pi8=8.d0*pi
      k_infinum=k_inf*k_fermi
      q_infinum=k_infinum
      r_infinum=500.d0
      n_q_heg=n_k_heg
      n_p_heg=4*n_k_heg
      n_x_heg=4*n_k_heg
      n_r_heg=4*n_k_heg
      alf_k_heg=0.05d0
      q_small_heg=0.004d0*k_fermi
      call heg_mesh
      call memory_heg
c --------- Coulomb interaction -------------------------
      v_q_heg=0.d0      
      do iq=1,n_q_heg
        if(dim_heg==2) v_q_heg(iq)=4.d0*pi/q_m_heg(iq)
        if(dim_heg==3) v_q_heg(iq)=8.d0*pi/q_m_heg(iq)**2
      enddo
      v_q_heg(0)=v_q_heg(1)
      if(iter_psi+iter_bsp/=0) then
        vol0=pi2**3/amega/nqdiv_c
        qs=(3.d0*vol0/4.d0/pi)**(1.d0/3.d0)
        v_q_heg_c(1)=32.d0*pi*pi*qs/vol0
        do i=2,nstar_c
          v_q_heg_c(i)=8.d0*pi/(pi2a*qabs_star_c(i))**2
        enddo
        call from_q_to_rr_heg_red(nstar_c,v_q_heg_c,v_r_heg_c)
      endif
c ----------------------------------------------------------------
      chem_pot=0.d0
      it0=iter_dft+iter_hf
      call timel('*** HEG SCF started ****************')
      iter=1
      ubi_0='dft'
      ro_spin=ro_el_gas/2
      sigx_k_heg=0.d0
      sigc_k_omega_heg=(0.d0,0.d0)
      sigc_k_omega_as=(0.d0,0.d0)
      do while (iter<=iter_max_scf)
        if(maswrk) write(iun,1010)iter
        call ubi_sumus(iter)        
        last=.false.
        if(ubi=='dft'.and.iter==iter_dft) then
          last=.true.
        else if(ubi==' hf'.and.iter==iter_dft+iter_hf) then
          last=.true.
        else if(ubi==' gw'.and.iter==iter_dft+iter_hf+iter_gw) then
          last=.true.
        else if(ubi==' qp') then
          if(iter==iter_dft+iter_hf+iter_gw+iter_qp) last=.true.
        else if(ubi=='psi') then
          if(iter==iter_dft+iter_hf+iter_gw+iter_qp+iter_psi) 
     &        last=.true.
        else if(ubi=='bsp') then
          if(iter==iter_dft+iter_hf+iter_gw+iter_qp+iter_psi+iter_bsp) 
     &        last=.true.
        endif
        call set_adm(adm_sig)
        call heg_b_split
c ----- Free Electron + External Magnetic Field spectrum ---------        
        do ispin=1,nspin
          w=b_split_heg*(-1)**(3-ispin)
          do i=0,n_k_heg
            e_star_free(i,ispin)=k_m_heg(i)**2+w
          enddo
        enddo
        if(ubi=='dft') then
          call heg_dft_spectrum
          if(iter==1) chem_pot=vxc_heg(1)+1.d0
        endif
c ---------------------------------------------------------------------
        if(ubi==' qp'.and.ubi/=ubi_0) e_star_qp=e_star_x
        if(ubi/='dft'.and.ubi/=' hf') then
c ----------- P + W evaluation ----------------------------------------
          call polarizability_heg
          call output_heg1_w_nu
          call output_heg_w_tau
          call output_heg1_e_nu
          call output_heg1_resp_nu(1)
c	    call output_heg_e1_rax
c          call local_field_factors_1
c          call heg_vrt_q0_nu0(w_q_tau_heg)
c ----- Correlated part of Self Energy -----------------------------
          call heg1_sig_c(w_q_tau_heg)
          if(ubi=='psi'.or.ubi=='bsp') then
c -------------- Sigma-Vertex evaluation ------------------------------
            if(ubi=='psi') call heg_vertex_sig(psi_sig,iter_sigma_gwg)
            if(ubi=='bsp') call heg_vertex_sig(psi2_sig,iter_sigma_gwg)
            allocate(pl(nc_line))
            allocate(spl(nc_line-1,4))
            hh=(qabs_line_c(2)-qabs_line_c(1))*pi2a
            do ispin=1,nspin
              do i_tau=0,n_tau
                pl=sig_k_tau_vertex_heg(i_tau,:,ispin)
                call spline_hmg(hh,pl,spl(1,1),spl(1,2),spl(1,3),
     &                          spl(1,4),nc_line-1,0,0.d0,0.d0)
                do iq=0,n_k_heg
                  q=k_m_heg(iq)
                  if(q>hh*(nc_line-1)) then
                    sig_k_tau_vertex_rad(iq,i_tau,ispin)=0.d0
                  else
                    i=q/hh
                    if(i*hh<q) i=i+1
                    if(i==0) i=1
                    q=q-hh*(i-1)
                    sig_k_tau_vertex_rad(iq,i_tau,ispin)=
     &                                   spl(i,1)+q*(spl(i,2)+q*(
     &                                   spl(i,3)+q*spl(i,4)))
                  endif
                enddo
                do i=0,n_k_heg
                  sigc_k_tau_heg(i,i_tau,ispin)=
     &              sigc_k_tau_heg(i,i_tau,ispin)
     &             +sig_k_tau_vertex_rad(i,i_tau,ispin)
                enddo
              enddo
            enddo
            deallocate(pl,spl)
            do ispin=1,nspin
              call tau_to_omega_spl(
     &                          sig_k_omega_vertex_heg(0,0,ispin),
     &	                      sig_k_tau_vertex_rad(0,0,ispin),
     &	                            dsigc_k_tau_heg(0,1,ispin),
     &                                n_k_heg+1,dsig_tau_to_omega)
            enddo
          endif
          if(ubi==' qp') then
            do ispin=1,nspin
              call tau_to_omega_spl(sigc_qp_k_omega_heg(0,0,ispin),
     &	                          sigc_k_tau_heg(0,0,ispin),
     &	                          dsigc_k_tau_heg(0,1,ispin),
     &                              n_k_heg+1,dsig_tau_to_omega)
            enddo
          endif
c ------------- Mixing the Self Energy -------------------------         
          if(iter>iter_dft+iter_hf+1) sigc_k_tau_heg=
     &	    adm_sig*sigc_k_tau_heg+(1.d0-adm_sig)*sigc_k_tau_heg_old
          do ispin=1,nspin
            call tau_to_omega_spl(sigc_k_omega_heg(0,0,ispin),
     &	                        sigc_k_tau_heg(0,0,ispin),
     &	                        dsigc_k_tau_heg(0,1,ispin),
     &                            n_k_heg+1,dsig_tau_to_omega)
          enddo
        endif
c ---------- New X-spectrum ----------------------------------------
        if(ubi/='dft') then
c ------- Effective exchange -----------------------------------
          call heg1_sig_x
          call output_heg1_sig_x
c ----------- Mix Sigma_X -------------------------------	
          if(iter>iter_dft+1) sigx_k_heg=adm_sig*sigx_k_heg
     &	                              +(1.d0-adm_sig)*sigx_k_heg_old	
          e_star_x=e_star_free+sigx_k_heg
        endif
        if(ubi/='dft'.and.ubi/=' hf') then
          call output_heg1_sig_c_tau
          call output_heg1_sig_omega
        endif
c ------------- New Chemical Potential ------------------------       
        call search_mu_0(chem_pot,1.d0,chem_pot)
        if(maswrk) write(iun,*)' Chemical Potential: ',chem_pot
        if(ubi=='dft') chem_pot_lda=chem_pot
c -------- New Exchange Green's Function -------------------------------
        gx_k_tau_heg=0.d0
        gx_k_omega_heg=(0.d0,0.d0)
        do ispin=1,nspin
         do i=0,n_k_heg
            de=e_star_x(i,ispin)-chem_pot
            do i_tau=0,n_tau
              gx_k_tau_heg(i,i_tau,ispin)=green0(de,tau_mesh(i_tau))
            enddo
            dgx_k_tau_heg(i,1,ispin)=-de*gx_k_tau_heg(i,0,ispin)
            dgx_k_tau_heg(i,2,ispin)=-de*gx_k_tau_heg(i,n_tau/2,ispin)
            dgx_k_tau_heg(i,3,ispin)=-de*gx_k_tau_heg(i,n_tau,ispin)
            do i_omega=0,n_omega
              gx_k_omega_heg(i,i_omega,ispin)=(1.d0,0.d0)
     &	                           /dcmplx(-de,w_omega(i_omega))
            enddo
          enddo
        enddo
        g_k_tau_heg=gx_k_tau_heg
        g_k_omega_heg=gx_k_omega_heg
        dg_k_tau_heg=dgx_k_tau_heg
c -------- New Correlated part of Green's Function ------------
        if(ubi==' qp') then
          do ispin=1,nspin
            do i=0,n_k_heg
              de=e_star_qp(i,ispin)-chem_pot
              do i_omega=0,n_omega
                g_k_omega_heg(i,i_omega,ispin)=(1.d0,0.d0)
     &	                              /dcmplx(-de,w_omega(i_omega))
                gc_k_omega_heg(i,i_omega,ispin)=
     &	            g_k_omega_heg(i,i_omega,ispin)
     &	           -gx_k_omega_heg(i,i_omega,ispin)
              enddo
              do i_tau=0,n_tau
                g_k_tau_heg(i,i_tau,ispin)=green0(de,tau_mesh(i_tau))
                gc_k_tau_heg(i,i_tau,ispin)=g_k_tau_heg(i,i_tau,ispin)
     &	                           -gx_k_tau_heg(i,i_tau,ispin)
              enddo
            enddo   !! over i
          enddo
        else if(ubi/='dft'.and.ubi/=' hf') then
          call bands_qp_heg1(chem_pot)
          do ispin=1,nspin
            do i=0,n_k_heg
              do i_omega=0,n_omega
                gxx=gx_k_omega_heg(i,i_omega,ispin)
                gxs=gxx*sigc_k_omega_heg(i,i_omega,ispin)
                gc_k_omega_heg(i,i_omega,ispin)=gxs*gxx
     &	                                     /((1.d0,0.d0)-gxs)
              enddo
            enddo   !! over i
            call omega_to_tau(n_k_heg+1,
     &                        gc_k_omega_heg(0,0,ispin),
     &                        gc_k_tau_heg(0,0,ispin),
     &                        dgc_k_tau_heg(0,1,ispin))
          enddo
          g_k_omega_heg=g_k_omega_heg+gc_k_omega_heg
          g_k_tau_heg=g_k_tau_heg+gc_k_tau_heg
          dg_k_tau_heg=dg_k_tau_heg+dgc_k_tau_heg
        endif
        call heg_splines
c ----------- G_R_Tau ----------------------------------------
        call heg1_g_r_tau
        call output_heg1_g_tau
        ro_spin(1)=-g_r_tau_heg(0,n_tau,1)
        ro_spin(2)=-g_r_tau_heg(0,n_tau,nspin)
        if(nspin==2.and.maswrk) then
          write(iun,*)' Density of Spin Magnetization is ',
     &                ro_spin(2)-ro_spin(1)
          write(iun,*)' Spin Magnetization is ',
     &                (ro_spin(2)-ro_spin(1))*amega
          if(b_extval>1.d-8) write(iun,*)
     &     ' Uniform Spin Susceptibility ((Bm)^2/Ry/a.u.v.) is ',
     &                (ro_spin(2)-ro_spin(1))/b_extval
        endif
c        call output_heg_g_tau
c -----------------------------------------------------------------
        if(ubi/='dft'.and.ubi/=' hf') then
          call output_heg1_z_ren
        endif
c --------- Internal energy calculation --------------------------- 
        call etot_heg
c -----------------------------------------------------------------        
        if(last) then
          if(ubi=='dft') then
            call heg1_ra(rs_el_gas,0.d0)
c            call heg_ppm
c            call output_heg1_p_nu(p0_q_nu_heg,2)
          endif
          if(ubi=='dft'.or.ubi==' hf') then
            deltae=(emaxdos-emindos)/ndos
            do ispin=1,nspin
              do i_omega=0,ndos
                w=emindos+i_omega*deltae
                c1=dcmplx(w,e_small)
                do i=0,n_k_heg
                  g_rax_heg(i_omega,i,ispin)=(1.d0,0.d0)
     &		                        /(c1+chem_pot-e_star_x(i,ispin))
                enddo
              enddo
            enddo
          else if(ubi==' qp') then
            deltae=(emaxdos-emindos)/ndos
            do ispin=1,nspin
              do i_omega=0,ndos
                w=emindos+i_omega*deltae
                c1=dcmplx(w,e_small)
                do i=0,n_k_heg
                  g_rax_heg(i_omega,i,ispin)=(1.d0,0.d0)
     &		                        /(c1+chem_pot-e_star_qp(i,ispin))
                enddo
              enddo
            enddo
          endif
          if(ubi/='dft'.and.ubi/=' hf') then
            call g_rax_heg_corr
          endif
          call output_heg1_g_rax
          call output_sig_heg1_real
          call eff_mass_heg1
          call heg_susceptibility
        endif
        if(ubi/='dft') sigx_k_heg_old=sigx_k_heg
        if(ubi/='dft'.and.ubi/=' hf') sigc_k_tau_heg_old=sigc_k_tau_heg
        if(iter==iter_dft) chem_pot_lda=chem_pot
        iter=iter+1
        ubi_0=ubi
        if(iter>iter_max_scf) exit
      enddo   !!! over iter
c ---------------------------------------------------------------      
1010  format(/'                iteration ',i3/)
      call timel('**** HEG_RUN finished **************')
      end
