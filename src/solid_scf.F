      subroutine solid_scf
	  use atom_mod
	  use etot_mod
	  use manager_mod
	  use parallel_mod
	  use scalapack_mod
	  use solid_mod
	  use vertex_mod
      implicit none
      logical :: correlation,new_space
	  integer :: iexit1
      maxqmt=maxel*(maxel+1)/2
      call rostart(r_atom,ro_atom)
	  call init_vh
	  call writero
      call timel('***** rostart finished *************')
*  *****  the loop of self-consistency *****
      call timel('*******   SCF started **************')
	  iter=1
	  iexit1=0
      w_red_print=.false.
	  do while (iter<=iter_max_scf)
	    call ubi_sumus(iter)
        key_hybrid=0
        if(ubi==' hf'.and.hybrid_type/=0) key_hybrid=1
	    if(ubi/='dft'.and.ubi/=' hf') then
	      if(lsclp.and.iter==iter_dft+iter_hf+1) then
c -------- Prepare for SCALAPACK using in W evaluations ------------
#ifdef SCLP
           call sclp_pw_prep
#endif
          endif
	      if(restart_begin>restart_end) then
	        iter=restart_end
	        call ubi_sumus(iter)
	        call restart_solid(1)
            goto 1
	      endif
	      if(restart_begin>iter) then
	        iter=restart_begin
	        call ubi_sumus(iter)
	        call restart_solid(1)
	      endif
	    endif
        if(maswrk) write(iun,1010)iter
	    if(ubi=='dft') ubi_0=ubi
	    if(nspin/=nspin_1) nspmag=max(magn,nspin)
c ---------------------------------------------------------------------     
        we_need_pol=.false.
        if(ubi==' qp') we_need_pol=.true.
        if(ubi==' gw') we_need_pol=.true.
        if(ubi=='psi') we_need_pol=.true.
        if(ubi=='bsp') we_need_pol=.true.
	    if(iter/=restart_begin) call set_g00
c ---------------------------------------------------------------------
        if(.not.iheg) call smultro(ro,rointr,v_intr_h,vatmt)
        call vcoul(ro,vatmt,v_mt_h,v_intr_h)
        call add_vxc     
	    call e_v_xch(1,vxch) !! for old Density
        call core_all
	    call output_ro_v
        if(itest/=0.and..not.iheg) then
	      text=' Matching : Coulon Potential!'
          call matching(0,v_mt_h,v_intr_h)
        endif
	    if(ubi=='dft') then
	      call dft_loop
        else if(ubi==' gw'.and.mbpt_gw/=0) then
          call mbpt_loop
        else
	      call gw_loop
	      call metal_gw
	    endif
	    call metal_mf
c ------------------------------------------------------------------
        call e_v_xch(2,exch_dft)
	    correlation=.false.
        if(ubi==' qp'.or.ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') 
     &    correlation=.true.
        call ROfull(correlation,ro,rointr)
	    call coulener(ro_vh_old)
        if(.not.iheg) call smultro(ro,rointr,v_intr_h,vatmt)
        call vcoul(ro,vatmt,v_mt_h,v_intr_h)
	    call coulener(ro_vh_new)
*     admix density for the next iteration
        if(iter>=1.and.ubi=='dft') call mixro1
        if(hybrid_type==2.or.hybrid_type==4) then
          if(hybrid_type==2) alpha_hbr=1.10507d0*gbar_param-1.00778d0
          if(hybrid_type==4) alpha_hbr=0.130711d0*gbar_param**4
     &                                +0.121983d0
          if(maswrk) write(iun,*)' Alpha_hbr =',alpha_hbr
        endif
	    call writero
	    iexit1=0
	    if(ubi=='dft'.and.iter==iter_dft) iexit1=1
	    if(ubi==' hf'.and.iter==iter_dft+iter_hf) iexit1=1
	    if(ubi==' gw'.and.iter==iter_dft+iter_hf+iter_gw) iexit1=1
	    if(ubi==' qp'.and.iter==iter_dft+iter_hf+iter_gw
     &                                          +iter_qp) iexit1=1
	    if(ubi=='psi'.and.iter==iter_dft+iter_hf+iter_gw
     &                                  +iter_qp+iter_psi) iexit1=1
	    if(ubi=='bsp'.and.iter==iter_dft+iter_hf+iter_gw
     &                         +iter_qp+iter_psi+iter_bsp) iexit1=1
	    if(iter==iter_dft) then
          we_need_v=.false.
          if(iter_max_scf/=iter_dft) we_need_v=.true.
          if(chi_cmp(1:1)=='1') we_need_v=.true.
          if(opt_mode(1:1)/='0'.and.opt_loc_fld) we_need_v=.true.
          if(iter_max_scf/=iter_dft.or.chi_cmp/='000'.or.
     &       (opt_mode(1:1)/='0'.and.opt_loc_fld)) then
	        call prod_bas
	        if(l_double_freq) call prod_bas_red
            if(vrt_x_appr(2:2)=='2') call trans_pb_pw_long
	      endif
	    endif
	    if(iexit1==1) then
c ---- Low energy bands for vertex corrections --------------
          if(nbndf_bnd_0(2)/=0) call bands_vertex
	      new_space=.false.
	      if(ubi=='dft') new_space=.true.
	      if(new_space) then
            if(nbndf_bnd_0(2)/=0.and.l_double_freq) then
              call bands_vertex1
c ---- Integrals for PW-susceptibilities ----------------------------
              call psi_psi_pw_red
            endif
c --------------------------------------------------------------	 
	      endif
	    endif
        if(iter==iter_dft) then
          allocate(e_lda(nbndf,npnt,nspin))
          e_lda=e_bnd
          chem_pot_lda=chem_pot
        endif
1	    call properties(iexit1)
	    if(iter==restart_end.and.restart_begin<=restart_end) then
	      call restart_solid(-1)
          if(lsclp) then
c -------- Deallocate SCALAPACK arrays used in W evaluations -------
#ifdef SCLP
            call sclp_pw_fin
#endif
          endif
	      exit
	    endif
        if(ubi/='dft'.and.ubi/=' hf') then
          if(lsclp.and.iter==iter_max_scf) then
c -------- Deallocate SCALAPACK arrays used in W evaluations -------
#ifdef SCLP
            call sclp_pw_fin
#endif
          endif
        endif
	    ubi_0=ubi
 	    iter=iter+1
      enddo   !!! over iter
1010  format(/'                iteration ',i3/)
      end
