      subroutine susceptibility_dft
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical :: cht
      integer :: ind_ch0(2),ind_ix0(3)
c -------- LDA - based responses ----------------------------------
      if(chi_cmp/='000') then
        cht=.false.
        if(vrt_x_appr(3:3)=='3') cht=.true.
        call resp_lda_prepare(chi_cmp,cht)
        if(nspin==1) then
          if(vrt_x_appr(1:1)=='1') then
            if(chi_cmp(1:1)=='1') call resp_lda_uniform('___00',1,1)
            if(chi_cmp(3:3)=='1') call resp_lda_uniform('___ZZ',1,1)
          endif
          if(vrt_x_appr(2:2)=='2') then
            if(chi_cmp(1:1)=='1') call resp_lda_static('___00',1,1)
          endif
          if(vrt_x_appr(3:3)=='3') then
            if(chi_cmp(1:1)=='1') call resp_lda('___00',1,1,1,1)
            if(chi_cmp(3:3)=='1') call resp_lda('___ZZ',1,1,1,1)
          else if(vrt_x_appr(3:3)=='1') then
            if(chi_cmp(1:1)=='1') call resp_lda_q1('___00',1,1)
            if(chi_cmp(3:3)=='1') call resp0_lda_q1_pw('___ZZ')
            if(chi_cmp(3:3)=='1') call resp_lda_q1('___ZZ',1,1)
          endif
        else if(nspin==2) then
          if(chi_cmp(1:1)=='1'.or.chi_cmp(3:3)=='1') then
            ind_ch0(1)=ind_chi(1)
            ind_ch0(2)=ind_chi(2)
            ind_ix0(1)=1
            ind_ix0(2)=2
            ind_ix0(3)=3
            if(.not.allocated(p_head_0z)) then
              allocate(p_head_0z(3,3,ndim3_nu))
              allocate(p_wing_0z(3,n_opt_pb,ndim3_nu))
              allocate(aa01_0z(n_opt_pb))
              call p_head_wings_dft(p_head_0z,p_wing_0z,aa00_0z,aa01_0z,
     &                              '0Z')
            endif
            if(vrt_x_appr(3:3)=='3') then
              call resp_lda('00_ZZ',2,ind_ch0,3,ind_ix0)
            else if(vrt_x_appr(3:3)=='1') then
              call resp_lda_q1('00_ZZ',2,ind_ch0,3,ind_ix0)
            endif
          endif
          if(chi_cmp(2:2)=='1') then
            ind_ch0(1)=ind_chi(3)
            ind_ch0(2)=ind_chi(4)
            ind_ix0(1)=n_ixc_0
            if(vrt_x_appr(3:3)=='3') then
              call resp_lda('___XY',2,ind_ch0,1,ind_ix0)
            else if(vrt_x_appr(3:3)=='1') then
              call resp0_lda_q1_pw('___XY')
              call resp_lda_q1('___XY',1,1)
            endif
          endif
        endif
        deallocate(theta_int)
        if(chi_basis=='PB') deallocate(thet_me)
        if(chi_basis=='PW') deallocate(theta_mt)
        if(cht) deallocate(resp_lda_tau_mm,resp_lda_tau_mi,
     &                     resp_lda_tau_ii)
      endif
      call timel('**** SUSCEPTIBILITY_DFT finished ***')
      end
