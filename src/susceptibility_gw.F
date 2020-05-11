      subroutine susceptibility_gw
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      if(vrt_x_appr/='000') then
        allocate(w_red_q(n_pbtot_red,n_pbtot_red,ndim3_nu,npnt))
        call w_c_from_b
! ------- Uniform susceptibilities (00 and ZZ) ------------------
        if(vrt_x_appr(1:1)=='1') then
          if(chi_cmp(1:1)/='0'.or.chi_cmp(3:3)/='0')
     &                                  call response_uniform_00_zz
        endif
! ------- Static (nu=0) susceptibilities along Q_LINE --------------
c ------- Uses full PB_red basis for evaluation of Xi, E, E^-1 -----
        if(vrt_x_appr(2:2)=='2') then
          if(chi_cmp(1:1)/='0'.or.chi_cmp(3:3)/='0')
     &                                  call response_static_00_zz
        endif
! ------- Dynamic (nu>=0) susceptibilities at one Q-point ----------
! ------- Uses only PW (G=G'=0) basis ------------------------------
        if(vrt_x_appr(3:3)=='1') then
          if(chi_cmp(1:1)/='0') call response_dynamic_pw(1)
          if(chi_cmp(2:2)/='0') call response_dynamic_pw(2)
          if(chi_cmp(3:3)/='0') call response_dynamic_pw(3)
        endif
! ------- Dynamic (nu>=0) susceptibilities at one Q-point ----------
c ------- Uses full PB_red basis for evaluation of Xi, E, E^-1 -----
        if(vrt_x_appr(3:3)=='1') then
          if(chi_cmp(1:1)/='0') call response_full
        endif
        deallocate(w_red_q)
      endif
      call timel('**** SUSCEPTIBILITY_GW finished ****')
      end
