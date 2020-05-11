      subroutine set_solid
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
      integer :: isort,jsort
c     limlb - number of different L,m or L,i,mj
      if(irel.ne.2) then
        limlb1=(maxb+2)**2
      else if(irel.eq.2) then
        limlb1=2*(maxb+2)**2
      endif
      limlb=nrel*(maxb+1)**2
      limlps=(maxs+1)**2
      limlb_red=nrel*(maxb_red+1)**2
      limlpb=(maxpb+1)**2  !! for product basis
c -----------------------------------------------------------------
      allocate(iob(natom))
      allocate(iob2(natom))
      allocate(iobs(natom))
      allocate(iob20(nsort))
      allocate(amass(nsort))
      allocate(lfun(nsort))
      allocate(lfunm(nsort))
      allocate(ncor(nsort))
      allocate(ncorm(nsort))
      allocate(lc_max(nsort))
      allocate(io_lem(natom))
      allocate(io_lem0(nsort))
      allocate(io_cor(natom))
      allocate(nlctot(nsort))
      allocate(indbas(maxntle,limlb,natom))
      allocate(indfun0(2,maxntle,nrel*maxb+1,nsort))
      call afteread
      allocate(indbasa(2,maxntle,limlb,nsort))
      allocate(lm_isz(maxel,nsort))
      allocate(lf_isz(maxel,nsort))
      allocate(ie_isz(maxel,nsort))
      allocate(indcor(ncormax,nspin_0,nsort))
      indcor=0
      allocate(ind_wf(maxlfun,nsort))
      ind_wf=0
      allocate(ntle_cor(nrel*maxb+1,nsort))
      allocate(indcor0(7,nrel*maxb+1,nsort))
      allocate(nod_cor(ncormax,nsort))
      allocate(li_cor(ncormax,nsort))
      allocate(ie_val(maxlfun,nsort))
      allocate(li_val(maxlfun,nsort))
      allocate(in_val(maxlfun,nsort))
      allocate(lm_cor(ncormmax,nsort))
      allocate(lf_cor(ncormmax,nsort))
      allocate(indcore(8,nrel*limlb,nsort))
      allocate(e_core(ncormax,nspin_0,nsort))
      allocate(psi_nre(ncormax,nspin_0,nsort))
      allocate(r_nre_core(ncormax,nspin_0,nsort))
      allocate(cor_norm(ncormax,nspin_0,nsort))
      call flapw_1
      if(mesh_k_c) call orbitals_red
      maxplw=800000
      nfgtl=5
      maxplwro=100000
      numvol=50
! ---- Here is for AFM irel<2 symmetrization  ---------------------
      allocate(iafm(nsort))
      iafm=0
      if(irel<2.and.nspin==2) then
        do isort=1,nsort
          if(magn_shift(isort)<1.d-3) cycle
          do jsort=isort+1,nsort
            if(abs(z(isort)-z(jsort))>1.d-3) cycle
            if(abs(magn_shift(isort)+magn_shift(jsort))<1.d-3) then
              iafm(isort)=jsort
              iafm(jsort)=isort
            endif
          enddo
        enddo
      endif
! ----------------------------------------------------------------
      mode_pb=1
      q0_long=(/0.d0,0.d0,0.d0/)
      q1_long=(/1.d0,0.d0,0.d0/)
      q_long_factor=6
      p_reducible=.false.
      mode_bts=0
      if(iter_gw/=0.and.w_sc_gw=='bts') mode_bts=1
      if(iter_qp/=0.and.w_sc_qp=='bts') mode_bts=1
      end
