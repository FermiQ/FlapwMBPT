      subroutine etot_heg
      use etot_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i_omega,i,i_nu,ispin
      real*8 :: t_kin,e_x,e_c_w,e_xc,etot,e1,etot1,e,e_c_pw,e_vrt_sg,
     &          e_vrt_pw,e_vrt_sg_om,de,fermi_dirac,sig_qp,e_c_qp,k2,
     &          ffd,e_vrt_nu,free_energy_x,wg,phic,ts
      t_kin=0.d0
      do ispin=1,nspin
        do i=0,n_k_heg
          k2=k_m_heg(i)**2
          ffd=g_k_tau_heg(i,n_tau,ispin)
          t_kin=t_kin-k2*k2*ffd*k_w_heg(i)
        enddo  !! over i
      enddo
      t_kin=t_kin/pi/pi/ro_el_gas/nspin
      e_x=0.d0
      if(ubi=='dft') then
        e_x=exc_heg(1)/2
        exch_dft=2.d0*e_x
      else
        d_free_x=0.d0
        do ispin=1,nspin
          do i=0,n_k_heg
            wg=k_m_heg(i)**2*k_w_heg(i)
            ffd=sigx_k_heg(i,ispin)*g_k_tau_heg(i,n_tau,ispin)
            e_x=e_x-wg*ffd
            ffd=free_energy_x(e_star_x(i,ispin)-chem_pot)
            d_free_x=d_free_x+wg*ffd
          enddo  !! over i
        enddo
        d_free_x=d_free_x/pi/pi/ro_el_gas/nspin
        e_x=e_x/4.d0/pi/pi/ro_el_gas/nspin
      endif
      e_x=2.d0*e_x    !! spin
c -------------------------------------------------------------
      e_vrt_sg=0.d0
      e_vrt_sg_om=0.d0
      e_vrt_pw=0.d0      
      e_c_w=0.d0
      e_c_pw=0.d0
      e_c_qp=0.d0
	if(ubi=='psi'.or.ubi=='bsp') then
        e_vertex_pw_tau0=0.d0
        e_vrt_nu=0.d0
      endif
      if(ubi/='dft'.and.ubi/=' hf') then
	  call heg_free_energy_c
	  grand_pot_x=d_free_x-e_x
	  grand_pot=grand_pot_x
	  phi_functional_c=0.d0
        do i=0,n_q_heg
          wg=q_m_heg(i)**2*q_w_heg(i)/pi/pi/ro_el_gas
	    do i_nu=0,n_nu
	      call heg_phi_rpa_c(p_q_nu_heg(i,i_nu),v_q_heg(i),phic)
            phi_functional_c=phi_functional_c
     &                      +et_sum_nu(i_nu)*phic*wg
	    enddo
	  enddo
        do ispin=1,nspin
          do i=0,n_k_heg
            k2=k_m_heg(i)**2
            if(ubi==' qp') then
              sig_qp=e_star_qp(i,ispin)-e_star_x(i,ispin)
              e_c_qp=e_c_qp-k2*k_w_heg(i)*g_k_tau_heg(i,n_tau,ispin)
     &                                   *sig_qp
            endif
            do i_omega=0,n_omega
              e=sigc_k_omega_heg(i,i_omega,ispin)
     &         *gx_k_omega_heg(i,i_omega,ispin)
              e1=sigc_k_omega_heg(i,i_omega,ispin)
     &          *gc_k_omega_heg(i,i_omega,ispin)
              e_c_w=e_c_w+(e*et_sum_sc_gx(i_omega)
     &                    +e1*et_sum_sc_gc(i_omega))*k2*k_w_heg(i)
              if(ubi=='psi'.or.ubi=='bsp') then
                e=sig_k_omega_vertex_heg(i,i_omega,ispin)
     &           *gx_k_omega_heg(i,i_omega,ispin)
                e1=sig_k_omega_vertex_heg(i,i_omega,ispin)
     &            *gc_k_omega_heg(i,i_omega,ispin)
                e_vrt_sg_om=e_vrt_sg_om+(e*et_sum_sc_gx(i_omega)
     &                                  +e1*et_sum_sc_gc(i_omega))
     &                             *k2*k_w_heg(i)
              endif
            enddo
          enddo
        enddo
        do i=0,n_q_heg
          k2=q_m_heg(i)**2
          do i_nu=0,n_nu
            e=p_q_nu_heg(i,i_nu)*w_q_nu_heg(i,i_nu)
            e_c_pw=e_c_pw-e*et_sum_nu(i_nu)*k2*q_w_heg(i)
          enddo
          if(ubi=='psi'.or.ubi=='bsp') then
            e=p_q_tau_vertex_heg(i,0)*v_q_heg(i)
            e_vertex_pw_tau0=e_vertex_pw_tau0-e*k2*q_w_heg(i)
            do i_nu=0,n_nu
              e=p_q_nu_vertex_rad(i,i_nu)*w_q_nu_heg(i,i_nu)
              e_vrt_nu=e_vrt_nu-e*et_sum_nu(i_nu)*k2*q_w_heg(i)
            enddo
          endif
        enddo
        if(ubi==' qp') then
          e_c_qp=e_c_qp/2.d0/pi/pi/ro_el_gas
        else
          e_c_w=e_c_w/pi/pi/ro_el_gas
          e_c_pw=e_c_pw/2.d0/pi/pi/ro_el_gas
        endif
        if(ubi=='psi'.or.ubi=='bsp') then
          e_vrt_sg_om=e_vrt_sg_om/pi/pi/ro_el_gas
          e_vertex_pw_tau0=e_vertex_pw_tau0/4.d0/pi/pi/ro_el_gas
          e_vrt_nu=e_vrt_nu/2.d0/pi/pi/ro_el_gas
          psi_functional_vrt=e_vrt_sg_om
          phi_functional_c=phi_functional_c+e_vertex_pw_tau0
        endif
	  grand_pot=grand_pot_x+phi_functional_c-2.d0*e_c_w+d_free_c
      endif
      if(ubi=='psi'.or.ubi=='bsp') then
        e_c_pw=e_c_pw+e_vertex_pw_tau0
        e_vrt_pw=e_vrt_nu+e_vertex_pw_tau0
      endif
      if(ubi/=' qp') e_xc=e_x+e_c_w
      if(ubi==' qp') e_xc=e_x+e_c_qp
      etot=t_kin+e_xc
      etot1=t_kin+e_x+e_c_pw
      if(ubi=='dft') then
c ---------- RPA(LDA) approach --------------------------------
        e_x_dft=0.d0
        do ispin=1,nspin
          do i=0,n_k_heg
            k2=k_m_heg(i)**2
            de=k2+vxc_heg(ispin)-chem_pot_lda
            e_x_dft=e_x_dft+k_w_heg(i)*k2*sigx_k_heg(i,ispin)
     &                                   *fermi_dirac(de)
          enddo
        enddo
        e_x_dft=e_x_dft/2.d0/pi/pi/ro_el_gas
        if(maswrk) then
          write(iun,*)' T_kin',t_kin
          write(iun,*)' E_X_DFT',e_x_dft
        endif
      endif
      if(maswrk) then
        write(iun,140) t_kin,e_x,e_vrt_sg,e_vrt_sg_om,e_vertex_pw_tau0,
     &                 e_vrt_pw,e_c_w,e_xc,e_c_pw,e_x+e_c_pw,ubi,etot,
     &                 etot1
        print *,ubi,': ETOT = ',etot,etot1
        write(iun,*)ubi,': ETOT = ',etot
        if(ubi/='dft') then
          write(iun,*)' D_Free_X ',d_free_x
          if(ubi/=' hf') then
	      ts=etot1-grand_pot-chem_pot
            write(iun,*)' D_Free_C ',d_free_c
            write(iun,*)' Phi_c_RPA ',phi_functional_c
            if(ubi=='psi'.or.ubi=='bsp') then
              write(iun,*)' Psi_c_VRT_Sig_G ',psi_functional_vrt
            endif
          endif
          write(iun,*)' Grand Potential ',grand_pot
          write(iun,*)' Free Energy F ',grand_pot+chem_pot
	    write(iun,*)' Entropy (Ry/K) ',ts/temperature
        endif
      endif
  140 format(/' T_KIN  = ',f15.8,' E_X = ',f15.8,
     &       /'Sig*G_Vertex_tau = ',f15.8,' Sig*G_Vertex_om = ',f15.8,
     &       /'P_VRT_tau0*V    = ',f15.8,' P*W_Vertex_full  = ',f15.8,
     &       /'SG_tau   E_C  = ',f15.8,'  E_XC_SG   = ',f15.8,
     .       /'PW_nu :  E_C  = ',f15.8,'  E_XC_PW   = ',f15.8,
     &       /a3,': ETOT (OM)  = ',f17.8,': ETOT_PW (OM)  = ',f17.8)
      end