      subroutine p_im_r(ind_tau,isp1,isp2,pmi)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer,intent(in) :: ind_tau,isp1,isp2
      real*8, intent(inout) :: pmi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                             ndim3_tn,ndim3_k(me_k+1))
      integer :: lm,l,ii,iii,indx,lget,icg1,icg2,icg,jatom,jsort,irl,
     &           km2,lm2,lf2,km1,lm1,lf1,i,ikm1,ikm2,ir,lf,ie,ind_ir,
     &           ispin1,ispin2
      real*8 :: cgg,dum
      complex*16 :: cggc,dumc
      real*8, allocatable :: gr(:,:,:,:),p_mi(:,:,:)
      complex*16, allocatable :: gc(:,:,:,:,:)
      ispin1=min(isp1,nspin_1)
      ispin2=min(isp2,nspin_1)
      if(irel<=1) then
        allocate(gr(nr_red,nfun,ndim3_kk(me_k+1),2))
        call g_rs_mi_mem(gr(1,1,1,1),dumc,1,isp1,ind_tau)
        call g_rs_mi_mem(gr(1,1,1,2),dumc,2,isp2,ind_tau)
      else if(irel==2) then
        allocate(gc(nr_red,nfun,ndim3_kk(me_k+1),nrel,2))
        call g_rs_mi_mem(dum,gc,1,1,ind_tau)
        call g_rs_mi_mem(dum,gc(1,1,1,1,2),2,1,ind_tau)
      endif
      if(ind_tau==1) call timel('** G_RS_MI finished ****************')
      allocate(p_mi(nr_red,n_pbmt,ndim3_kk(me_k+1)))
      p_mi=0.d0
      do ind_ir=1,ndim3_kk(me_k+1)
        ir=n3_mpi_kk(me_k+1)+ind_ir
        do jatom=1,natom
          jsort=is(jatom)
          do km1=1,lfunm(jsort)
            lm1=lm_isz(km1,jsort)
            lf1=lf_isz(km1,jsort)
            ikm1=io_lem(jatom)+km1-1
            do km2=1,lfunm(jsort)
              ikm2=io_lem(jatom)+km2-1
              lm2=lm_isz(km2,jsort)
              lf2=lf_isz(km2,jsort)
              if(irel<=1) then
                ii = max0(lm1,lm2)
                iii = min0(lm1,lm2)
                indx = (ii*(ii-1))/2 + iii
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1) - 1
                do icg = icg1, icg2
                  lm = jcg(icg)
                  l=lget(lm)
                  if(l>lmpb(jsort)) cycle
                  do ie=1,ntle_pb(l,jsort)
                    lf=ind_prod0(ie,l,jsort)
                    cgg=fi0(lf2,lf1,lf,jsort,ispin1,ispin2)*cg(icg)
                    i=iopb(jatom)+ind_prod(ie,lm,jsort)-1
                    p_mi(:,i,ind_ir)=p_mi(:,i,ind_ir)
     &                              -gr(:,ikm1,ind_ir,1)
     &                              *gr(:,ikm2,ind_ir,2)*cgg
                  enddo
                enddo
              else if(irel==2) then
                icg1=indxcgr(lm2,lm1,1,0)
                icg2=indxcgr(lm2,lm1,2,0)
                do icg = icg1, icg2
                  lm=jcgr(icg)
                  l=lget(lm)
                  if(l>lmpb(jsort)) cycle
                  do ie=1,ntle_pb(l,jsort)
                    lf=ind_prod0(ie,l,jsort)
                    cggc=fi0(lf2,lf1,lf,jsort,1,1)*cgr(icg)
                    i=iopb(jatom)+ind_prod(ie,lm,jsort)-1
                    do irl=1,2
                      p_mi(:,i,ind_ir)=p_mi(:,i,ind_ir)
     &                                -gc(:,ikm1,ind_ir,irl,1)
     &                                *conjg(gc(:,ikm2,ind_ir,irl,2))
     &                                *cggc
                    enddo
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      enddo  !! over ind_ir
      if(irel<=1) deallocate(gr)
      if(irel==2) deallocate(gc)
      if(ind_tau==1) call timel('** P_IM_R finished *****************')
      call p_im_q_from_r(p_mi,ind_tau,pmi)
      deallocate(p_mi)
      if(ind_tau==1) call timel('** P_IM_Q_FROM_R finished **********')
      end



      subroutine p_im_sclp(ind_tau,isp1,isp2,pwsclp)
      use atom_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      implicit none
      integer,intent(in) :: ind_tau,isp1,isp2
      real*8, intent(inout) :: pwsclp(ncmpl*nrow_pw_max,ncol_pw_max,
     &                                ndim3_tn,ndim3_k(me_k+1))
      integer :: lm,l,ii,iii,indx,lget,icg1,icg2,icg,jatom,jsort,irl,
     &           km2,lm2,lf2,km1,lm1,lf1,i,ikm1,ikm2,ir,lf,ie,ind_ir,
     &           ispin1,ispin2
      real*8 :: cgg,dum
      complex*16 :: cggc,dumc
      real*8, allocatable :: gr(:,:,:,:),p_mi(:,:,:)
      complex*16, allocatable :: gc(:,:,:,:,:)
      ispin1=min(isp1,nspin_1)
      ispin2=min(isp2,nspin_1)
      if(irel<=1) then
        allocate(gr(nr_red,nfun,ndim3_kk(me_k+1),2))
        call g_rs_mi_mem(gr(1,1,1,1),dumc,1,isp1,ind_tau)
        call g_rs_mi_mem(gr(1,1,1,2),dumc,2,isp2,ind_tau)
      else if(irel==2) then
        allocate(gc(nr_red,nfun,ndim3_kk(me_k+1),nrel,2))
        call g_rs_mi_mem(dum,gc,1,1,ind_tau)
        call g_rs_mi_mem(dum,gc(1,1,1,1,2),2,1,ind_tau)
      endif
      if(ind_tau==1) call timel('** G_RS_MI finished ****************')
      allocate(p_mi(nr_red,n_pbmt,ndim3_kk(me_k+1)))
      p_mi=0.d0
      do ind_ir=1,ndim3_kk(me_k+1)
        ir=n3_mpi_kk(me_k+1)+ind_ir
        do jatom=1,natom
          jsort=is(jatom)
          do km1=1,lfunm(jsort)
            lm1=lm_isz(km1,jsort)
            lf1=lf_isz(km1,jsort)
            ikm1=io_lem(jatom)+km1-1
            do km2=1,lfunm(jsort)
              ikm2=io_lem(jatom)+km2-1
              lm2=lm_isz(km2,jsort)
              lf2=lf_isz(km2,jsort)
              if(irel<=1) then
                ii = max0(lm1,lm2)
                iii = min0(lm1,lm2)
                indx = (ii*(ii-1))/2 + iii
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1) - 1
                do icg = icg1, icg2
                  lm = jcg(icg)
                  l=lget(lm)
                  if(l>lmpb(jsort)) cycle
                  do ie=1,ntle_pb(l,jsort)
                    lf=ind_prod0(ie,l,jsort)
                    cgg=fi0(lf2,lf1,lf,jsort,ispin1,ispin2)*cg(icg)
                    i=iopb(jatom)+ind_prod(ie,lm,jsort)-1
                    p_mi(:,i,ind_ir)=p_mi(:,i,ind_ir)
     &                              -gr(:,ikm1,ind_ir,1)
     &                              *gr(:,ikm2,ind_ir,2)*cgg
                  enddo
                enddo
              else if(irel==2) then
                icg1=indxcgr(lm2,lm1,1,0)
                icg2=indxcgr(lm2,lm1,2,0)
                do icg = icg1, icg2
                  lm=jcgr(icg)
                  l=lget(lm)
                  if(l>lmpb(jsort)) cycle
                  do ie=1,ntle_pb(l,jsort)
                    lf=ind_prod0(ie,l,jsort)
                    cggc=fi0(lf2,lf1,lf,jsort,1,1)*cgr(icg)
                    i=iopb(jatom)+ind_prod(ie,lm,jsort)-1
                    do irl=1,2
                      p_mi(:,i,ind_ir)=p_mi(:,i,ind_ir)
     &                                -gc(:,ikm1,ind_ir,irl,1)
     &                                *conjg(gc(:,ikm2,ind_ir,irl,2))
     &                                *cggc
                    enddo
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      enddo  !! over ind_ir
      if(irel<=1) deallocate(gr)
      if(irel==2) deallocate(gc)
      if(ind_tau==1) call timel('** P_IM_R finished *****************')
      call p_im_q_sclp(p_mi,ind_tau,pwsclp)
      deallocate(p_mi)
      if(ind_tau==1) call timel('** P_IM_Q_FROM_R finished **********')
      end
