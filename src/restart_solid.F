      subroutine restart_solid(key)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: key
      logical :: w_read
      integer :: i_len,n
#ifdef MPI
      integer :: isum
#endif
c ------------------------------------------------------------
      i_len=len_trim(tmpfile)
      if(key.gt.0) then
        if(maswrk) then
          open(3,file=tmpfile(1:i_len)//'_solid_rst',form='unformatted')
          read(3)chem_pot,ubi_0,w_read
          read(3)e_core,pcor,qcor
          read(3)n_bnd,e_bnd
          read(3)g_loc_0
          read(3)gfun,p_f,pd_f,pd2_f,q_f,qd_f,qd2_f
          if(irel>=1) read(3)gfund
          read(3)eny,ffsmt
          if(mode_bts==1) read(3)fxc_00_bts
          close(3)
        endif
        if(goparr) then
#ifdef MPI
          call brdcst(msgdbl,chem_pot,8,master,MPI_COMM_WORLD)
          call brdcstc(ubi_0,3,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,e_core,8*ncormax*nspin*nsort,master,
     &	            MPI_COMM_WORLD)
          call brdcst(msgdbl,pcor,8*maxmtcor,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,qcor,8*maxmtcor,master,MPI_COMM_WORLD)
          call brdcst(msgint,w_read,4*npnt*nspin,master,MPI_COMM_WORLD)
          call brdcst(msgint,n_bnd,4*npnt*nspin,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,e_bnd,8*nbndf*npnt*nspin,master,
     &	            MPI_COMM_WORLD)
          call brdcst(msgdbl,g_loc_0,16*maxel*maxel*natom*nspin,master,
     &	            MPI_COMM_WORLD)
          call brdcst(msgdbl,gfun,8*maxwf*nspin,master,MPI_COMM_WORLD)
          if(irel>=1) call brdcst(msgdbl,gfund,8*maxwf*nspin,master,
     &	                        MPI_COMM_WORLD)
          isum=(maxnrad+1)*maxntle*(nrel*maxb+1)*nsort*nspin
          call brdcst(msgdbl,p_f,8*isum,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,pd_f,8*isum,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,pd2_f,8*isum,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,q_f,8*isum,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,qd_f,8*isum,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,qd2_f,8*isum,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,eny,8*maxntle*(nrel*maxb+1)*nsort*nspin,
     &	            master,MPI_COMM_WORLD)
          call brdcst(msgdbl,ffsmt,8*(maxntle+1)**2*(nrel*maxb+1)
     &		          *nsort*nspin,master,MPI_COMM_WORLD)
          if(mode_bts==1) call brdcst(msgdbl,fxc_00_bts,8*n_phead_dir,
     &                                master,MPI_COMM_WORLD)
#endif     
        endif
! -------------- G_FULL ----------------------------------------
        if(me_b==0) call write_g_corr(1,1,g_full)
        if(nproc_b/=1) then
          n=16*nbndf**2*ndim3_tau*ndim3_k(me_k+1)*nspin
          call brdcst(msgdbl,g_full,n,0,comm_b)
        endif
! -------------- SIG_C_TAU_OLD ----------------------------------------
        if(me_b==0) call write_g_corr(1,2,sig_c_tau_old)
        if(nproc_b/=1) then
          n=16*nbndf**2*ndim3_tau*ndim3_k(me_k+1)*nspin
          call brdcst(msgdbl,sig_c_tau_old,n,0,comm_b)
        endif
! -------------- V_H_BB_OLD ----------------------------------------
        call write_sig_x(1,1,v_h_bb_old)
! -------------- SIGX_SOLID_K_OLD -----------------------------------
        call write_sig_x(1,2,sigx_solid_k_old)
! -------------- EV_BND --------------------------------------------
        if(me_t_b==0) call write_a_bnd(1)
        if(nproc_t*nproc_b/=1) then
          n=16*nbasmpw*nbndf*ndim3_k(me_k+1)*nspin
          call brdcst(msgdbl,ev_bnd,n,0,comm_t_b)
        endif
! -------------- Z_BND ---------------------------------------------
        if(me_t_b==0) call write_z_bnd(1)
        if(nproc_t*nproc_b/=1) then
          n=16*nfun*nbndf*ndim3_k(me_k+1)*nspin
          call brdcst(msgdbl,z_bnd,n,0,comm_t_b)
        endif
! ------------ W_RED_Q --------------------------------------------
        if(w_read) then
          if(me_b==0) call write_w_red(1,w_nu_q_b)
          if(nproc_b/=1) then
            n=8*n_pbtot_red**2*ndim3_nu*ndim3_k(me_k+1)
            call brdcst(msgdbl,w_nu_q_b,n,0,comm_b)
          endif
          if(.not.allocated(w11_red))
     &      allocate(w11_red(n_opt_pb_red,n_opt_pb_red,ndim3_nu))
          if(me_k_b==0) call write_w11_red(1)
          if(nproc_k*nproc_b/=1) then
            n=16*n_opt_pb_red**2*ndim3_nu
            call brdcst(msgdbl,w11_red,n,0,comm_k_b)
          endif
        endif
        if(.not.allocated(ws_head)) allocate(ws_head(0:n_tau/2))
        if(.not.allocated(ws_head_nu)) allocate(ws_head_nu(0:n_nu))
        if(ncmpl==1) then
          if(.not.allocated(ws01_met_r)) allocate(ws01_met_r(n_opt_pb))
        else
          if(.not.allocated(ws01_met)) allocate(ws01_met(n_opt_pb))
        endif
        if(me==0) then
          call write_ws_head(1)
          if(iter>iter_dft+iter_hf+iter_gw+iter_qp+1)
     &       call write_vrt_scale(1)
        endif
#ifdef MPI
        if(nproc/=1) then
          call brdcst(msgdbl,ws_head,4*n_tau+8,master,MPI_COMM_WORLD)
          call brdcst(msgdbl,ws_head_nu,8*n_nu+8,master,MPI_COMM_WORLD)
          if(ncmpl==1) then
            call brdcst(msgdbl,ws01_met_r,8*n_opt_pb,master,
     &                  MPI_COMM_WORLD)
          else
            call brdcst(msgdbl,ws01_met,16*n_opt_pb,master,
     &                  MPI_COMM_WORLD)
          endif
          if(iter>iter_dft+iter_hf+iter_gw+iter_qp+1)
     &       call brdcst(msgdbl,vrt_scale,16*nbndf_bnd**2*npnt_c*nspin,
     &                   master,MPI_COMM_WORLD)
        endif
#endif
        if(ncmpl==1) then
          if(.not.allocated(ws_head_tens_r))
     &       allocate(ws_head_tens_r(3,3,ndim3_nu))
          if(.not.allocated(ws_wing_r)) allocate(ws_wing_r(3,n_opt_pb,
     &                                                     ndim3_nu))
        else
          if(.not.allocated(ws_head_tens)) allocate(ws_head_tens(3,3,
     &                                                        ndim3_nu))
          if(.not.allocated(ws_wing)) allocate(ws_wing(3,n_opt_pb,
     &                                                 ndim3_nu))
        endif
        if(me_k_b==0) call write_ws_head_par(1)
        if(nproc_k*nproc_b/=1) then
          if(ncmpl==1) then
            n=72*ndim3_nu
            call brdcst(msgdbl,ws_head_tens_r,n,0,comm_k_b)
            n=24*n_opt_pb*ndim3_nu
            call brdcst(msgdbl,ws_wing_r,n,0,comm_k_b)
          else
            n=144*ndim3_nu
            call brdcst(msgdbl,ws_head_tens,n,0,comm_k_b)
            n=48*n_opt_pb*ndim3_nu
            call brdcst(msgdbl,ws_wing,n,0,comm_k_b)
          endif
        endif
        if(mode_bts==1) call write_bts(1)
        pcor_old=pcor
        qcor_old=qcor
        gfun_old=gfun
        if(irel>=1) gfund_old=gfund
        call read_ro
        if(nbndf_bnd_0(2)/=0) call bands_vertex
        if(mesh_k_c) then
c ---- Low energy bands for vertex corrections --------------
          call bands_vertex1
c ---- Integrals for PW-susceptibilities ----------------------------
          call psi_psi_pw_red
        endif
        call timel('**** RESTART_SOLID: READING is OK **')
      else if(key.eq.-1) then
        if(maswrk) then
          open(3,file=tmpfile(1:i_len)//'_solid_rst',form='unformatted')
          write(3)chem_pot,ubi,w_red_print
          write(3)e_core,pcor,qcor
          write(3)n_bnd,e_bnd
          write(3)g_loc_0
          write(3)gfun,p_f,pd_f,pd2_f,q_f,qd_f,qd2_f
          if(irel>=1) write(3)gfund
          write(3)eny,ffsmt
          if(mode_bts==1) write(3)fxc_00_bts
          close(3)
        endif
! -------------- G_FULL ----------------------------------------
        if(me_b==0) call write_g_corr(0,1,g_full)
! -------------- SIG_C_TAU_OLD ----------------------------------------
        if(me_b==0) call write_g_corr(0,2,sig_c_tau_old)
! -------------- V_H_BB_OLD ----------------------------------------
        call write_sig_x(0,1,v_h_bb_old)
! -------------- SIGX_SOLID_K_OLD ----------------------------------
        call write_sig_x(0,2,sigx_solid_k_old)
! -------------- EV_BND --------------------------------------------
        if(me_t_b==0) call write_a_bnd(0)
! -------------- Z_BND ---------------------------------------------
        if(me_t_b==0) call write_z_bnd(0)
! ------------ W_RED_Q --------------------------------------------
        if(w_red_print) then
          if(me_b==0) call write_w_red(0,w_nu_q_b)
          if(me_k_b==0) call write_w11_red(0)
        endif
        if(me==0) then
          call write_ws_head(0)
          if(ubi=='psi'.or.ubi=='bsp') call write_vrt_scale(0)
        endif
        if(me_k_b==0) call write_ws_head_par(0)
        if(mode_bts==1) call write_bts(0)
        call write_ro
        call timel('**** RESTART_SOLID: WRITING is OK **')
      endif
      end
