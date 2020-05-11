      subroutine add_sigma_vertex_bands(ps)
c ------ Add correction to the full self energy ---------------
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ps
      integer :: ispin,ind_k,k,ind_tau,it,nbnd,nb,i,k1,k10,jb,ib,ka
      real*8 :: de,const
      real*8, allocatable :: sig(:,:,:,:),tmp2(:,:,:)
      complex*16, allocatable :: tmp1(:,:,:),tmp(:,:,:),a(:,:)
      e_vertex_sig_g_bnd=0.d0
      if(ps=='00000') return
      allocate(sig(nbndf,nbndf,2,ndim3_tau)) 
      allocate(tmp(nbndf_bnd,nbndf_bnd,2))
      allocate(tmp1(nbndf,nbndf,2))
      allocate(tmp2(nbndf,nbndf,2))
      allocate(a(nbndf_bnd,nbndf))
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          nbnd=n_bnd(k,ispin)
          sig=0.d0
          do ind_tau=1,ndim3_tau
            tmp1=(0.d0,0.d0)
            do i=1,num_baryc_ac(k)
              call zone1_number(pnt_baryc_ac(1,i,k),rb0_c,ndiv_c,k1)
              k1=index_k1_c(k1)   ! index in 0 BZ
              k10=i_kref_c(k1)
              ka=k_a_from_c(k10)
              nb=n_low_bnd(ka,ispin)
              call ferm_unpack_tau(tmp,
     &                             sig_bnd_vertex(1,1,1,ind_tau,k10,
     &                                            ispin),nb,
     &                             nbndf_bnd,nbndf_bnd,3)
              do it=1,2
                call zgemm('n','n',nb,nbnd,nb,(1.d0,0.d0),tmp(1,1,it),
     &                     nbndf_bnd,psps_ac(1,1,i,ind_k,ispin),
     &                     nbndf_bnd,(0.d0,0.d0),a,nbndf_bnd)
                call zgemm('c','n',nbnd,nbnd,nb,
     &	                   dcmplx(coef_baryc_ac(i,k),0.d0),
     &                     psps_ac(1,1,i,ind_k,ispin),nbndf_bnd,
     &	                   a,nbndf_bnd,(1.d0,0.d0),tmp1(1,1,it),nbndf)
              enddo
            enddo   !! over i
            tmp2=0.d0
            call ferm_pack_tau(tmp1,tmp2,nbnd,nbndf,nbndf)
            do jb=1,nbnd
              do ib=1,nbnd
                do it=1,2
                  sig_c_tau(ib,jb,it,ind_tau,ind_k,ispin)=
     &              sig_c_tau(ib,jb,it,ind_tau,ind_k,ispin)
     &             +tmp2(ib,jb,it)
                  sig(ib,jb,it,ind_tau)=tmp2(ib,jb,it)
                enddo
              enddo
            enddo
          enddo   !! over ind_tau
          call int_tau_sigc_gc(de,sig,g_full(1,1,1,1,ind_k,ispin),
     &                         nbnd,nbndf,nbndf)
          const=2.d0*wgt(k)/dfloat(nspin*nrel)
          e_vertex_sig_g_bnd=e_vertex_sig_g_bnd+const*de
        enddo   !! over ind_k
      enddo   !! over ispin
      if(nproc_k/=1) call dgop(e_vertex_sig_g_bnd,1,'  +',comm_k)
      if(nproc_t/=1) call dgop(e_vertex_sig_g_bnd,1,'  +',comm_t)
      deallocate(sig,tmp1,tmp2,tmp,a)
      end
