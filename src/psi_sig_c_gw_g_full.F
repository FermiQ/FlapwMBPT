      subroutine psi_sig_c_gw_g_full(psif)
c ------ Calculates Tr (Sig_C_gw * G_full) on fine k-mesh ----------
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      real*8, intent(out) :: psif
      integer :: ispin,ind_k,k,n
      real*8 :: de,const
      psif=0.d0
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          call int_tau_sigc_gc(de,sig_c_tau(1,1,1,1,ind_k,ispin),
     &                         g_full(1,1,1,1,ind_k,ispin),n,nbndf,
     &                         nbndf)
          const=4.d0*wgt(k)/dfloat(nspin*nrel)
          psif=psif+const*de
        enddo   !! over ind_k
      enddo   !! over ispin
      if(nproc_k/=1) call dgop(psif,1,'  +',comm_k)
      if(nproc_t/=1) call dgop(psif,1,'  +',comm_t)
      end
