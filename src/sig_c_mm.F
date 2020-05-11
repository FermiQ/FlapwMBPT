      subroutine sig_c_mm(ind_tau)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau
      integer :: ind_k,ispin,k,nrr_maxs,isort,jsort,nat2,n2pbm,n2pbmt,n,
     &           nf2,nfun2,nat1,n1pbm,n1pbmt,nf1,nfun1,it
      integer, allocatable :: r0_pairs(:,:,:),ig_pairs(:,:,:),
     &                        nrr_reds(:,:),nd_nrr(:),n0_nrr(:)
      real*8, allocatable :: w_rs(:,:,:)
      complex*16, allocatable :: s_c(:,:,:),pvt(:,:,:,:,:),tmp(:,:,:)
      allocate(pvt(ndim_b_nfun(me_b+1),nbndf,2,ndim3_k(me_k+1),nspin))
      allocate(nd_nrr(nproc_k))
      allocate(n0_nrr(nproc_k))
      pvt=(0.d0,0.d0)
      do jsort=1,nsort
        nat2=nhsort(jsort)
        n2pbm=n_pbmt0(jsort)
        n2pbmt=n2pbm*nat2
        nf2=lfunm(jsort)
        nfun2=nf2*nat2
        do isort=1,nsort
          nat1=nhsort(isort)
          n1pbm=n_pbmt0(isort)
          n1pbmt=n1pbm*nat1
          nf1=lfunm(isort)
          nfun1=nf1*nat1
          allocate(r0_pairs(nqdiv,nat1,nat2))
          allocate(ig_pairs(nqdiv,nat1,nat2))
          allocate(nrr_reds(3,nqdiv*nat1*nat2))
          call get_rr_sort(isort,jsort,nat1,nat2,r0_pairs,ig_pairs,
     &                     nrr_reds,nrr_maxs)
          call size_shift_par(nrr_maxs,nproc_k,nd_nrr,n0_nrr)
          allocate(w_rs(n1pbm,n2pbm,nd_nrr(me_k+1)))
          call w_mm_sort(ind_tau,w_rs,n1pbm,n2pbm,nd_nrr(me_k+1),
     &                   n0_nrr(me_k+1),n1pbmt,n2pbmt,nat1,nat2,isort,
     &                   jsort,nrr_reds)
          do ispin=1,nspin
            do it=1,2
              allocate(s_c(nfun1,nfun2,ndim3_kk(me_k+1)))
              call sigc_mm_r_sort(ispin,it,ind_tau,s_c,w_rs,n1pbm,n2pbm,
     &                            nd_nrr(me_k+1),n0_nrr(me_k+1),nfun1,
     &                            nfun2,isort,jsort,nf1,nf2,nrr_maxs,
     &                            nat1,nat2,nrr_reds,r0_pairs,ig_pairs)
              call sigc_mm_k_from_r_sort(ispin,it,s_c,nfun1,nfun2,
     &                                   pvt(1,1,1,1,ispin),isort,jsort)
              deallocate(s_c)
            enddo
          enddo
          deallocate(r0_pairs,ig_pairs,nrr_reds,w_rs)
        enddo  !! over isort
      enddo   !! over jsort
      deallocate(nd_nrr,n0_nrr)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          allocate(tmp(n,n,2))
          do it=1,2
            call zgemm('c','n',n,n,ndim_b_nfun(me_b+1),(1.d0,0.d0),
     &                 z_bnd(nmpi_b_nfun(me_b+1)+1,1,ind_k,ispin),nfun,
     &                 pvt(1,1,it,ind_k,ispin),ndim_b_nfun(me_b+1),
     &                 (0.d0,0.d0),tmp(1,1,it),n)
          enddo
          if(nproc_b/=1) call dgop(tmp,4*n*n,'  +',comm_b)
          call ferm_pack_tau(tmp,sig_c_tau(1,1,1,ind_tau,ind_k,ispin),n,
     &                       n,nbndf)
          deallocate(tmp)
        enddo
      enddo
      deallocate(pvt)
      end
