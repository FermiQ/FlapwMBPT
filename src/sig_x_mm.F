      subroutine sig_x_mm
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: ind_k,ispin,k,nrr_maxs,isort,jsort,nat2,n2pbm,n2pbmt,n,
     &           nf2,nfun2,nat1,n1pbm,n1pbmt,nf1,nfun1
      integer, allocatable :: r0_pairs(:,:,:),ig_pairs(:,:,:),
     &                        nrr_reds(:,:),nd_nrr(:),n0_nrr(:)
      real*8, allocatable :: v_rs_mm(:,:,:)
      complex*16, allocatable :: s_c(:,:,:),pvt(:,:,:,:),tmp(:,:)
      if(ndim_tb_nfun(me_t_b+1)/=0)
     & allocate(pvt(ndim_tb_nfun(me_t_b+1),nbndf,ndim3_k(me_k+1),nspin))
      allocate(nd_nrr(nproc_k))
      allocate(n0_nrr(nproc_k))
      if(ndim_tb_nfun(me_t_b+1)/=0) pvt=(0.d0,0.d0)
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
          if(maswrk) write(iun,*)' ISORT JSORT NRR_MAXS = ',isort,jsort,
     &                           nrr_maxs
          call size_shift_par(nrr_maxs,nproc_k,nd_nrr,n0_nrr)
          allocate(v_rs_mm(n1pbm,n2pbm,nd_nrr(me_k+1)))
          call v_mm_from_q_sort(v_rs_mm,n1pbm,n2pbm,nd_nrr(me_k+1),
     &                          n0_nrr(me_k+1),n1pbmt,n2pbmt,nat1,nat2,
     &                          isort,jsort,nrr_reds)
          allocate(pvv(nfun1,nfun2,ndim3_k(me_k+1)))
          do ispin=1,nspin
            pvv=(0.d0,0.d0)
            allocate(s_c(nfun1,nfun2,ndim3_kk(me_k+1)))
            call sigx_mm_r_sort(ispin,s_c,v_rs_mm,n1pbm,n2pbm,
     &                          nd_nrr(me_k+1),n0_nrr(me_k+1),nfun1,
     &                          nfun2,isort,jsort,nf1,nf2,nrr_maxs,nat1,
     &                          nat2,nrr_reds,r0_pairs,ig_pairs)
            call sigx_mm_k_from_r_sort(s_c,nfun1,nfun2)
            deallocate(s_c)
            if(ndim_tb_nfun(me_t_b+1)/=0) then
              do ind_k=1,ndim3_k(me_k+1)
                k=n3_mpi_k(me_k+1)+ind_k
                call sigma_mt_mt_fin_sort(ind_k,k,ispin,isort,jsort,
     &                                    nfun1,nfun2,
     &                                    pvt(1,1,ind_k,ispin))
              enddo
            endif
          enddo
          deallocate(r0_pairs,ig_pairs,nrr_reds,pvv,v_rs_mm)
        enddo  !! over isort
      enddo   !! over jsort
      deallocate(nd_nrr,n0_nrr)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          allocate(tmp(n,n))
          tmp=(0.d0,0.d0)
          if(ndim_tb_nfun(me_t_b+1)/=0) then
            call zgemm('c','n',n,n,ndim_tb_nfun(me_t_b+1),(1.d0,0.d0),
     &                 z_bnd(nmpi_tb_nfun(me_t_b+1)+1,1,ind_k,ispin),
     &                 nfun,pvt(1,1,ind_k,ispin),ndim_tb_nfun(me_t_b+1),
     &                 (0.d0,0.d0),tmp,n)
          endif
          if(nproc_t*nproc_b/=1) call dgop(tmp,2*n*n,'  +',comm_t_b)
          call pack_hrm_tb(tmp,sigx_solid_k(1,ind_k,ispin),n,n,1)
          deallocate(tmp)
        enddo
      enddo
      if(ndim_tb_nfun(me_t_b+1)/=0) deallocate(pvt)
      end
