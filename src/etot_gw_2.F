      subroutine etot_gw_2
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      logical :: last,grnd
      integer :: isort
      real*8 :: gvh,t_kin_x,t_kin_c,e_coul,t_kin,e_x,e_c,e_xc,etot,
     &	        etot_pw,ts,mu_n,free_x_val,ts_new,ts0,ncr
      ncr=0.d0
      do isort=1,nsort
        ncr=ncr+nhsort(isort)*zcor(isort)
      enddo
      last=.false.
      if(ubi==' hf'.and.iter==iter_dft+iter_hf) last=.true.
      if(ubi==' gw'.and.iter==iter_dft+iter_hf+iter_gw) last=.true.
      if(ubi=='psi'.and.iter==iter_dft+iter_hf+iter_gw+iter_qp+iter_psi)
     &                                                  last=.true.
      if(ubi=='bsp'.and.iter==iter_dft+iter_hf+iter_gw+iter_qp+iter_psi
     &                                                        +iter_bsp)
     &                                                  last=.true.
      grnd=.false.
      if(ubi/='dft'.and.key_hybrid==0) then
        if(psi_fncl_use==2.or.(psi_fncl_use==1.and.last)) grnd=.true.
        if(ubi==' gw'.and.w_sc_gw=='non') grnd=.false.
        if(ubi=='psi'.and.w_sc_psi=='non') grnd=.false.
        if(ubi=='bsp'.and.w_sc_bsp=='non') grnd=.false.
      endif
      if(ubi=='dft'.or.ubi_0=='dft'.or.key_hybrid==1) then
        gx_vxc_val=vxch-gx_vxc_cor
      endif
      if(nproc_k/=1) then
        call DGOP(gx_vh_val,1,'  +',comm_k)
        call DGOP(ebnd_x,1,'  +',comm_k)
        call DGOP(gx_sx_val,1,'  +',comm_k)
        call DGOP(d_free_x,1,'  +',comm_k)
        if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
          call DGOP(gc_vh_val,1,'  +',comm_k)
          call DGOP(ebnd_c,1,'  +',comm_k)
          call DGOP(gc_sx_val,1,'  +',comm_k)
          call DGOP(gc_sc_val,1,'  +',comm_k)
          call DGOP(gx_sc_val,1,'  +',comm_k)
          call DGOP(gc_sc_val_w,1,'  +',comm_k)
          call DGOP(gx_sc_val_w,1,'  +',comm_k)
          call DGOP(gc_vxc_val,1,'  +',comm_k)
          call DGOP(d_free_c,1,'  +',comm_k)
        endif
      endif
      if(nproc_t/=1) then
        if(ubi_0=='dft'.or.key_hybrid==1) 
     &	  call brdcst(msgdbl,gc_vxc_val,8,0,comm_t)
        if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
          call brdcst(msgdbl,ebnd_c,8,0,comm_t)
          call DGOP(gc_sc_val,1,'  +',comm_t)
          call DGOP(gx_sc_val,1,'  +',comm_t)
          call DGOP(gc_sc_val_w,1,'  +',comm_t)
          call DGOP(gx_sc_val_w,1,'  +',comm_t)
        endif
      endif
      if(maswrk)
     &write(iun,139) z_vnucl,ro_vh_old,ro_vh_new,ecor,ebnd_x,ebnd_c,
     &               gx_vxc_cor/2,gx_vxc_val/2,gc_vxc_val/2,
     &               gx_vh_cor/2,gx_vh_val/2,gc_vh_val/2,
     &               gx_sx_cor/2,gx_sx_val/2,gc_sx_val/2,
     &               gx_sc_val/2,gc_sc_val/2,gx_sc_val_w/2,
     &               gc_sc_val_w/2
  139 format(/' Z_VNUCL=',f17.8,' ROVHOLD=',f17.8,' ROVHNEW=',f17.8,
     &      /' ECOR=   ',f16.8,'   EBND_X= ',f15.8,'  EBND_C = ',f15.8,
     &     /'GX_VXC_COR=',f15.8,' GX_VXC_VAL=',f14.8,' GC_VXC_V=',f14.8,
     &      /'GX_VH_COR=',f16.8,' GX_VH_VAL=',f15.8,' GC_VH_VAL=',f15.8,
     &      /'GX_SX_COR=',f15.8,' GX_SX_VAL=',f15.8,' GC_SX_VAL=',f15.8,
     &      /' GX_SC_VAL = ',f15.8,' GC_SC_VAL = ',f15.8,
     &      /' GX_SC_VAL_W=',f15.8,' GC_SC_VAL_W=',f15.8)
      mu_n=chem_pot*(ncr+nelec)
      t_kin_c=0.d0
      if(ubi=='dft') then
        e_coul=z_vnucl+0.5d0*ro_vh_new
      else
        gvh=gx_vh_cor+gx_vh_val+gc_vh_val
        e_coul=z_vnucl+0.5d0*gvh
      endif
      t_kin_x=ecor+ebnd_x
      if(ubi=='dft') then
        t_kin_x=t_kin_x-ro_vh_old-vxch
      else 
        t_kin_x=t_kin_x-(gx_vh_cor+gx_vh_val)
        t_kin_x=t_kin_x-gx_sx_cor-gx_sx_val
        if(key_hybrid==1) t_kin_x=t_kin_x-vxch
      endif
      if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
        t_kin_c=ebnd_c-gc_vh_val-gc_sx_val
      endif
c -------------- Convolutions SIGMA*GREEN -----------------------
      e_x=0.d0
      e_c=0.d0
      if(grnd) then
        free_x_val=d_free_x
        ts0=0.d0
        d_free_x=d_free_x+ecor-chem_pot*ncr
      endif
      if(ubi=='dft') then
        e_xc=exch_dft
      else
        e_x=0.5d0*(gx_sx_cor+gx_sx_val+gc_sx_val)
        if(key_hybrid==1) e_x=e_x+exch_dft
        if(grnd) then
          grand_pot_x=d_free_x+e_coul-gvh-e_x
          grand_pot=grand_pot_x
        endif
        if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
          e_c=0.5d0*(gx_sc_val+gc_sc_val)
          if(grnd) then
            grand_pot=grand_pot_x+phi_functional_c-2.d0*e_c+d_free_c
            psi_functional_c=psi_c_1_0
            if(maswrk) write(iun,*)'E_C_SOLID_PW',e_c_solid_pw
            if(ubi=='psi'.or.ubi=='bsp') then
              if(maswrk) then
                write(iun,*)'PSI_P_vrt_Ws',psi_p_vrt_ws
                write(iun,*)'E_vertex_PW_Tau0',e_vertex_pw_tau0
              endif
              psi_functional_c=psi_c_1_0+psi_functional_vrt
              grand_pot=grand_pot+psi_functional_vrt+psi_p_vrt_ws
            endif
          endif
        endif
        e_xc=e_x+e_c
      endif
      t_kin=t_kin_x+t_kin_c
      etot=t_kin+e_coul+e_xc+e_b_ext
      if(maswrk.and.(ubi=='psi'.or.ubi=='bsp')) then
        write(iun,*)' E_vertex_sig_g_FULL_BND ',e_vertex_sig_g_bnd
        if(grnd) then
          write(iun,'(a20,e12.5,a6,e12.5)')' Psi_c_1_0:   Fine:',
     &                             psi_c_1_0
          write(iun,'(a20,e12.5,a6,e12.5)')' Psi_c_1_1: Coarse:',
     &                             psi_c_1_1,' Fine:',psi_c_1_1_fine
          write(iun,'(a24,e12.5)')' Psi_c_1_1 (PW): Coarse:',
     &                          -0.5d0*(psi_p_vrt_ws-e_vertex_pw_tau0)
          write(iun,*)' Psi_Functional_vrt',psi_functional_vrt
          write(iun,*)' Psi_Functional_c',psi_functional_c
        endif
      endif
      if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
        if(grnd) then
          if(ubi=='psi') then
            if(w_sc_psi=='scf') e_c_solid_pw=e_c_solid_pw
     &                                    +e_vertex_pw_tau0
            if(w_sc_psi/='scf'.and.ubi/=ubi_0) e_c_solid_pw=
     &                               e_c_solid_pw+e_vertex_pw_tau0
          endif
          if(ubi=='bsp') then
            if(w_sc_bsp=='scf') e_c_solid_pw=e_c_solid_pw
     &                                    +e_vertex_pw_tau0
            if(w_sc_bsp/='scf'.and.ubi/=ubi_0) e_c_solid_pw=
     &                               e_c_solid_pw+e_vertex_pw_tau0
          endif
          if(maswrk) then
            etot_pw=t_kin+e_coul+e_x+e_c_solid_pw+e_b_ext
            write(iun,140) e_coul,t_kin,e_x,e_c,e_xc,e_c_solid_pw,etot,
     &                     etot_pw
            if(grnd) then
              print '(a3,1x,a20,3f15.8)',ubi,': E(SG), E(PW), F = ',
     &                                 etot,etot_pw,grand_pot+mu_n
            else
              print '(a3,1x,a17,2f19.12)',ubi,': E(SG), E(PW) = ',
     &                              etot,etot_pw
            endif
          endif
        else if(maswrk) then
          write(iun,141) e_coul,t_kin,e_x,e_c,e_xc,etot
          print '(a3,1x,a10,f19.12)',ubi,': E(SG) = ',etot
        endif
      else if(maswrk) then
        write(iun,141) e_coul,t_kin,e_x,e_c,e_xc,etot
        if(grnd) print *,ubi,': ETOT F = ',etot,grand_pot+mu_n
        if(.not.grnd) print *,ubi,': ETOT = ',etot
      endif
      if(grnd) then
        if(ubi==' hf') ts=etot-grand_pot-mu_n
        if(ubi/=' hf') then
          ts=etot-grand_pot-mu_n
          ts0=ebnd_c-d_free_c+3.d0*e_c-phi_functional_c
        endif
        ts_new=ebnd_x-free_x_val+ts0-chem_pot*nelec
      endif
      if(maswrk) then
        if(grnd) then
          write(iun,*)' D_Free_X, D_Free_C ',d_free_x,d_free_c
          write(iun,*)' Phi_c_RPA ',phi_functional_c
          write(iun,*)' Grand Potential ',grand_pot
          write(iun,*)' Free Energy F ',grand_pot+mu_n
          write(iun,*)' Entropy (Ry/K) ',ts/temperature
          if(ubi/='dft') write(iun,*)' Entropy (Ry/K)N',
     &                   ts_new/temperature
        endif
      endif
      call timel('***** ETOT_GW_2 finished ***********')
  140 format(/' E_COUL  = ',f17.8,' T_KIN  = ',f15.8
     .       /'GS : E_X = ',f15.8,'   E_C  = ',f15.8,'  E_XC   = ',f15.8
     .       /'PW : E_C  = ',f15.8
     .       /' ETOT (GS)  = ',f17.8
     .       /' ETOT (PW)  = ',f17.8)
  141 format(/' E_COUL  = ',f17.8,' T_KIN  = ',f15.8
     .       /'GS : E_X = ',f15.8,'   E_C  = ',f15.8,'  E_XC   = ',f15.8
     .       /' ETOT (GS)  = ',f17.8)
      end
