      subroutine psi_sig_g_vrt(sig,psif)
c ------ Calculates Tr (Sig*Gfull) on coarse k-mesh ----------------
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      complex*16, intent(in) :: sig(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                              ndim_k_red(me_k+1),nspin)
      real*8, intent(out) :: psif
      integer :: ispin,ind_k,k,ind_tau,it,n,kc,ibnd,jbnd,ib,jb
      real*8 :: de,const
      real*8, allocatable :: g_tau(:,:,:,:,:),sig_tau(:,:,:,:)
      psif=0.d0
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c))
      allocate(sig_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau))
      do ispin=1,nspin
        g_tau=0.d0
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          kc=k_c_from_a(k)
          if(kc==0) cycle
          n=n_low_bnd(k,ispin)
          do ind_tau=1,ndim3_tau
            do it=1,2
              do jbnd=1,n
                jb=ind_bands_bnd(jbnd,k,ispin)
                do ibnd=1,n
                  ib=ind_bands_bnd(ibnd,k,ispin)
                  g_tau(ibnd,jbnd,it,ind_tau,kc)=
     &                g_full(ib,jb,it,ind_tau,ind_k,ispin)
                enddo
              enddo
            enddo
          enddo   !! over ind_tau
        enddo  ! over ind_k
        if(nproc_k/=1) call dgop(g_tau,2*nbndf_bnd**2*npnt_c*ndim3_tau,
     &                           '  +',comm_k)
        do ind_k=1,ndim_k_red(me_k+1)
          k=n_mpi_k_red(me_k+1)+ind_k
          n=n_low_bnd(k,ispin)
          sig_tau=(0.d0,0.d0)
          do ind_tau=1,ndim3_tau
            call ferm_pack_tau(sig(1,1,1,ind_tau,ind_k,ispin),
     &                         sig_tau(1,1,1,ind_tau),n,nbndf_bnd,
     &                         nbndf_bnd)
          enddo   !! over ind_tau
          call int_tau_sigc_gc(de,sig_tau,g_tau(1,1,1,1,k),n,nbndf_bnd,
     &                         nbndf_bnd)
          const=4.d0*wgt_c(k)/dfloat(nspin*nrel)
          psif=psif+const*de
        enddo   !! over ind_k
      enddo   !! over ispin
      if(nproc_k/=1) call dgop(psif,1,'  +',comm_k)
      if(nproc_t/=1) call dgop(psif,1,'  +',comm_t)
      deallocate(g_tau,sig_tau)
      end
