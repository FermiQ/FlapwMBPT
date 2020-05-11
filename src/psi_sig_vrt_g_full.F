      subroutine psi_sig_vrt_g_full(sig,psif)
c ------ Calculates Tr (Sig_vrt * G_full) on fine k-mesh -----------
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
      integer :: ispin,ind_k,k,ind_tau,it,n,i,k1,k10,
     &           ka,nb
      real*8 :: de,const
      real*8, allocatable :: sigf(:,:,:,:,:)
      complex*16, allocatable :: tmp1(:,:,:),a(:,:),sg(:,:,:,:)
      psif=0.d0
      allocate(sigf(nbndf,nbndf,2,ndim3_tau,ndim3_k(me_k+1)))
      allocate(tmp1(nbndf,nbndf,2))
      allocate(a(nbndf_bnd,nbndf))
      allocate(sg(nbndf_bnd,nbndf_bnd,2,npnt_c))
      do ispin=1,nspin
        sigf=0.d0
        do ind_tau=1,ndim3_tau
          sg=(0.d0,0.d0)
          do it=1,2
            do ind_k=1,ndim_k_red(me_k+1)
              k=n_mpi_k_red(me_k+1)+ind_k
              n=n_low_bnd(k,ispin)
              sg(1:n,1:n,it,k)=sig(1:n,1:n,it,ind_tau,ind_k,ispin)
            enddo
          enddo
          if(nproc_k/=1) call dgop(sg,4*nbndf_bnd**2*npnt_c,'  +',
     &                             comm_k)
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            n=n_bnd(k,ispin)
            tmp1=(0.d0,0.d0)
            do i=1,num_baryc_ac(k)
              call zone1_number(pnt_baryc_ac(1,i,k),rb0_c,ndiv_c,k1)
              k1=index_k1_c(k1)   ! index in 0 BZ
              k10=i_kref_c(k1)
              ka=k_a_from_c(k10)
              nb=n_low_bnd(ka,ispin)
              do it=1,2
                call zgemm('n','n',nb,n,nb,(1.d0,0.d0),sg(1,1,it,k10),
     &                     nbndf_bnd,psps_ac(1,1,i,ind_k,ispin),
     &                     nbndf_bnd,(0.d0,0.d0),a,nbndf_bnd)
                call zgemm('c','n',n,n,nb,
     &                     dcmplx(coef_baryc_ac(i,k),0.d0),
     &                     psps_ac(1,1,i,ind_k,ispin),nbndf_bnd,
     &                     a,nbndf_bnd,(1.d0,0.d0),tmp1(1,1,it),nbndf)
              enddo
            enddo   !! over i
            call ferm_pack_tau(tmp1,sigf(1,1,1,ind_tau,ind_k),n,nbndf,
     &                         nbndf)
          enddo   !! over ind_k
        enddo   !! over ind_tau
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          call int_tau_sigc_gc(de,sigf(1,1,1,1,ind_k),
     &                         g_full(1,1,1,1,ind_k,ispin),n,nbndf,
     &                         nbndf)
          const=4.d0*wgt(k)/dfloat(nspin*nrel)
          psif=psif+const*de
        enddo   !! over ind_k
      enddo   !! over ispin
      if(nproc_k/=1) call dgop(psif,1,'  +',comm_k)
      if(nproc_t/=1) call dgop(psif,1,'  +',comm_t)
      deallocate(sigf,tmp1,a,sg)
      end
