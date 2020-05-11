      subroutine prod_bas
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: limax,ndim0,k
      real*8 :: pi2
      integer, allocatable :: nl(:,:),ind_l(:,:,:)
      pi2=pi+pi
      limax=nrel*maxb+1
c	  key=0 - Run to find the dimensions
      allocate(nl(limax,limax))
      allocate(ind_l(limax,limax,maxpb+1))
      call prod_b1(ndim0,nl,ind_l,limax)
      if(maswrk) write(iun,*)' Max dimension in full PB search is ',
     &                       ndim0
      call prod_b2(ndim0,0,nl,ind_l,limax)
      if(allocated(ind_pb)) deallocate(ind_pb)
      allocate(ind_pb(maxlfpb,nsort))
      if(allocated(lval_pb)) deallocate(lval_pb)
      allocate(lval_pb(maxlfpb,nsort))
      if(allocated(fun_pb)) deallocate(fun_pb)
      allocate(fun_pb(maxpbr))
      if(allocated(dfun_pb)) deallocate(dfun_pb)
      allocate(dfun_pb(maxpbr))
      if(allocated(lm_pbmt)) deallocate(lm_pbmt)
      allocate(lm_pbmt(n_pbmtm,nsort))
      if(allocated(lf_pb)) deallocate(lf_pb)
      allocate(lf_pb(n_pbmtm,nsort))
      if(allocated(ind_prod0)) deallocate(ind_prod0)
      allocate(ind_prod0(max_e_pb,0:maxpb,nsort))
      if(allocated(ind_prod)) deallocate(ind_prod)
      allocate(ind_prod(max_e_pb,limlpb,nsort))
c	  key=1 - Run to fill the arrays
      call prod_b2(ndim0,1,nl,ind_l,limax)
      if(allocated(fi0)) deallocate(fi0)
      allocate(fi0(maxlfun,maxlfun,maxlfpb,nsort,nspin,nspin))
	  if(allocated(mi0)) deallocate(mi0)
      allocate(mi0(maxlfpb,nsort))
      if(hybrid_type>2) then
        if(allocated(mi0_b)) deallocate(mi0_b)
        allocate(mi0_b(maxlfpb,nsort))
      endif
	  if(allocated(mi4)) deallocate(mi4)
      allocate(mi4(max_e_pb,nsort))
      if(allocated(tmm_pb)) deallocate(tmm_pb)
      allocate(tmm_pb(max_e_pb,max_e_pb,0:maxpb,nsort))
      call get_fi0_gw
      call get_tmm(0,tmm_pb)
      if(hybrid_type>2) then
        allocate(tmm_pb_b(max_e_pb,max_e_pb,0:maxpb,nsort))
        call get_tmm(1,tmm_pb_b)
      endif
      n_pbtot=n_pbmt+nplw_gw
      call size_shift_par(n_pbmt*(n_pbmt+1)/2,nproc_b,
     &                    ndim_pbmt,n_mpi_pbmt)
      if(.not.allocated(ndim_pbtt)) allocate(ndim_pbtt(nproc_t))
      if(.not.allocated(n_mpi_pbtt)) allocate(n_mpi_pbtt(nproc_t))
      call size_shift_par(n_pbmt,nproc_t,ndim_pbtt,n_mpi_pbtt)
      if(.not.allocated(ndim_pbt)) allocate(ndim_pbt(nproc_t))
      if(.not.allocated(n_mpi_pbt)) allocate(n_mpi_pbt(nproc_t))
      call size_shift_par(n_pbtot,nproc_t,ndim_pbt,n_mpi_pbt)
      k=nproc_t*nproc_b
      if(.not.allocated(ndim_tb_n_pbmt)) allocate(ndim_tb_n_pbmt(k))
      if(.not.allocated(nmpi_tb_n_pbmt)) allocate(nmpi_tb_n_pbmt(k))
      call size_shift_par(n_pbmt,k,ndim_tb_n_pbmt,nmpi_tb_n_pbmt)
      if(.not.allocated(nd_tb_pbmt)) allocate(nd_tb_pbmt(k))
      if(.not.allocated(n0_tb_pbmt)) allocate(n0_tb_pbmt(k))
      call size_shift_par(n_pbmt*(n_pbmt+1)/2,k,nd_tb_pbmt,n0_tb_pbmt)
      if(.not.allocated(nd_b_pbmt)) allocate(nd_b_pbmt(nproc_b))
      if(.not.allocated(n0_b_pbmt)) allocate(n0_b_pbmt(nproc_b))
      call size_shift_par(n_pbmt,nproc_b,nd_b_pbmt,n0_b_pbmt)
      if(mode_bts==1) then
        allocate(nd_tb_fxc(k))
        allocate(n0_tb_fxc(k))
        call size_shift_par(n_pbtot*(n_pbtot+1)/2,k,nd_tb_fxc,n0_tb_fxc)
        allocate(fxc_bts(ncmpl*nd_tb_fxc(me_t_b+1),ndim3_k(me_k+1)))
        allocate(fxc_00_bts(n_phead_dir))
        fxc_00_bts=0.d0
        fxc_bts=0.d0
      endif
      if(maswrk) then
        write(iun,*)' The size of PB is (k-point list) :'
        write(iun,'(14i5)')(n_pbmt+nplwgw(k),k=1,npnt)
        write(iun,'(a22,i6,a1,i3,a1,a5,i6,a7,i6)')
     &	               'Product basis: MT-part',n_pbmt,'(',nlb2,')',
     &                   ' Int:',nplw_gw,' Total:',n_pbtot
      endif
c ------------- Coulomb Interaction matrix -------------------------
      call trans_pb_pw
      if(iter_max_scf>iter_dft+iter_hf) then
        call bz_q0_num_prep
      endif
      if(we_need_v) then
        call av_q_prim
        if(allocated(v_mm_3)) deallocate(v_mm_3)
        allocate(v_mm_3(ncmpl*nd_tb_pbmt(me_t_b+1),ndim3_k(me_k+1)))
        if(allocated(v_mi_3)) deallocate(v_mi_3)
        allocate(v_mi_3(ncmpl*ndim_tb_n_pbmt(me_t_b+1),nplw_gw,
     &                  ndim3_k(me_k+1)))
        if(allocated(v_ii_3)) deallocate(v_ii_3)
        allocate(v_ii_3(nplw_gw,ndim3_k(me_k+1)))
        call v_coulomb(v_mm_3,v_mi_3,v_ii_3,0)
        if(iter_hf/=0.and.hybrid_type>2) then
          if(allocated(v_mm_3_sr)) deallocate(v_mm_3_sr)
          allocate(v_mm_3_sr(ncmpl*nd_tb_pbmt(me_t_b+1),
     &                       ndim3_k(me_k+1)))
          if(allocated(v_mi_3_sr)) deallocate(v_mi_3_sr)
          allocate(v_mi_3_sr(ncmpl*ndim_tb_n_pbmt(me_t_b+1),nplw_gw,
     &                       ndim3_k(me_k+1)))
          if(allocated(v_ii_3_sr)) deallocate(v_ii_3_sr)
          allocate(v_ii_3_sr(nplw_gw,ndim3_k(me_k+1)))
          call v_coulomb(v_mm_3_sr,v_mi_3_sr,v_ii_3_sr,1)
        endif
        if(ncmpl==1) call coulomb_ev_r
        if(ncmpl==2) call coulomb_ev
      endif
      deallocate(nl,ind_l)
      call timel('**** prod_bas finished *************')
      end
