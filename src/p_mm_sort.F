      subroutine p_mm_sort(pmm,ind_tau,isp1,isp2)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: isp1,isp2,ind_tau
      real*8, intent(inout) :: pmm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                             ndim3_k(me_k+1))
      integer :: nrr_maxs,isort,jsort,nat2,n2pbm,n2pbmt,nf2,nfun2,nat1,
     &           n1pbm,n1pbmt,nf1,nfun1
      integer, allocatable :: r0_pairs(:,:,:),ig_pairs(:,:,:),
     &                        nrr_reds(:,:),nd_nrr(:),n0_nrr(:)
      allocate(nd_nrr(nproc_k))
      allocate(n0_nrr(nproc_k))
      do jsort=1,nsort
        nat2=nhsort(jsort)
        n2pbm=n_pbmt0(jsort)
        n2pbmt=n2pbm*nat2
        nf2=lfunm(jsort)
        nfun2=nf2*nat2
        do isort=1,jsort
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
          call p_mm_mem_0(pmm,ind_tau,isp1,isp2,n1pbm,n2pbm,
     &                    nd_nrr(me_k+1),n0_nrr(me_k+1),nfun1,nfun2,
     &                    isort,jsort,nf1,nf2,nrr_maxs,nat1,nat2,
     &                    nrr_reds,r0_pairs,ig_pairs)
          deallocate(r0_pairs,ig_pairs,nrr_reds)
        enddo  !! over isort
      enddo   !! over jsort
      deallocate(nd_nrr,n0_nrr)
      if(ind_tau==1) call timel('** P_MM_SORT finished **************')
      end


      subroutine p_mm_sclp(pwsclp,ind_tau,isp1,isp2)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: isp1,isp2,ind_tau
      real*8, intent(inout) :: pwsclp(ncmpl*nrow_pw_max,ncol_pw_max,
     &                                ndim3_tn,ndim3_k(me_k+1))
      integer :: nrr_maxs,isort,jsort,nat2,n2pbm,n2pbmt,nf2,nfun2,nat1,
     &           n1pbm,n1pbmt,nf1,nfun1
      integer, allocatable :: r0_pairs(:,:,:),ig_pairs(:,:,:),
     &                        nrr_reds(:,:),nd_nrr(:),n0_nrr(:)
      allocate(nd_nrr(nproc_k))
      allocate(n0_nrr(nproc_k))
      do jsort=1,nsort
        nat2=nhsort(jsort)
        n2pbm=n_pbmt0(jsort)
        n2pbmt=n2pbm*nat2
        nf2=lfunm(jsort)
        nfun2=nf2*nat2
        do isort=1,jsort
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
          call p_mm_sclp_0(pwsclp,ind_tau,isp1,isp2,n1pbm,n2pbm,
     &                     nd_nrr(me_k+1),n0_nrr(me_k+1),nfun1,nfun2,
     &                     isort,jsort,nf1,nf2,nrr_maxs,nat1,nat2,
     &                     nrr_reds,r0_pairs,ig_pairs)
          deallocate(r0_pairs,ig_pairs,nrr_reds)
        enddo  !! over isort
      enddo   !! over jsort
      deallocate(nd_nrr,n0_nrr)
      if(ind_tau==1) call timel('** P_MM_SORT finished **************')
      end
