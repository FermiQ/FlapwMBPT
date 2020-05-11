      subroutine etot_gw_c(ispin,k,ind_k,sigc,ndim)
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ndim,k,ispin,ind_k
      real*8, intent(in) :: sigc(ndim,ndim,2,ndim3_omega)
      integer :: j
      real*8 :: const,gcom1,conv_hh_tb,tmp
      real*8, allocatable :: gb(:,:)
      allocate(gb(ndim,ndim))
      const=2.d0*wgt(k)/dfloat(nspin*nrel)
c -------------- G_beta --------------------------------------
      if(me_t==0) gb=0.5d0*(gc_tau(1:ndim,1:ndim,2,1)
     &                     -gc_tau(1:ndim,1:ndim,1,1))
      if(nproc_t/=1) call brdcst(msgdbl,gb,8*ndim**2,0,comm_t)
      if(ubi_0=='dft') then
        gc_vxc_val=gc_vxc_val
     &            -const*conv_hh_tb(gb,v_xc_bb(1,ind_k,ispin),ndim)
      endif
      gc_vh_val=gc_vh_val
     &         -const*conv_hh_tb(gb,v_h_bb(1,ind_k,ispin),ndim)
      tmp=-const*conv_hh_tb(gb,sigx_solid_k(1,ind_k,ispin),ndim)
      gc_sx_val=gc_sx_val+tmp
      do j=1,ndim
        ebnd_c=ebnd_c-const*gb(j,j)*e_bnd(j,k,ispin)
      enddo
      call free_energy_c0(ispin,ind_k,sigc,ndim)
c -------- GC_SC_VAL calculation --------------------------------------
      call int_tau_sigc_gc(gcom1,sig_c_tau(1,1,1,1,ind_k,ispin),
     &                     gc_tau,ndim,nbndf,nbndf)
      tmp=const*gcom1*2.d0
      gc_sc_val=gc_sc_val+tmp
c -------- GX_SC_VAL calculation --------------------------------------
      call int_tau_sigc_gx(gcom1,sig_c_tau(1,1,1,1,ind_k,ispin),
     &                     e_bnd(1,k,ispin),ndim,nbndf)
      tmp=const*gcom1*2.d0
      gx_sc_val=gx_sc_val+tmp
c -------- GC_SC_VAL_W calculation -------------------------------------
      call sum_maz_sig_gc(gcom1,sigc,gc_omega,ndim,nbndf,ndim)
      gc_sc_val_w=gc_sc_val_w+const*gcom1*2.d0
c -------- GX_SC_VAL_W calculation -------------------------------------
      call sum_maz_sig_gx(gcom1,sigc,e_bnd(1,k,ispin),ndim,ndim)
      gx_sc_val_w=gx_sc_val_w+const*gcom1*2.d0
      deallocate(gb)
      end
