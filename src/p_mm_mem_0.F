      subroutine p_mm_mem_0(pmm,ind_tau,isp1,isp2,n1pbm,n2pbm,
     &                      nd_nrr,n0_nrr,nfun1,nfun2,isort,jsort,nf1,
     &                      nf2,nrr_maxs,nat1,nat2,nrr_reds,r0_pairs,
     &                      ig_pairs)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: isp1,isp2,ind_tau,n1pbm,n2pbm,
     &                       nd_nrr,n0_nrr,nfun1,nfun2,isort,jsort,nf1,
     &                       nf2,nrr_maxs,nat1,nat2,
     &                       nrr_reds(3,nqdiv*nat1*nat2),
     &                       r0_pairs(nqdiv,nat1,nat2),
     &                       ig_pairs(nqdiv,nat1,nat2)
      real*8, intent(inout) :: pmm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                             ndim3_k(me_k+1))
      integer :: lm3,lm,l,ii,iii,indx,lget,icg1,icg2,icg,km3,lf3,
     &           jatom,km2,lm2,lf2,l1m1,l1,km1,lm1,lf1,iatom,
     &           km4,lm4,lf4,ie,lf,km,ir0,ind_ir0,ind_km3,ndim0,ni,
     &           nj,nipb,njpb,ispi1,ispi2
      real*8 :: cgg,eps,g2
      complex*16 :: g2c,cggc
      integer, allocatable :: ndim_km3(:),index_km3(:,:),nd_p2(:),
     &                        n0_p2(:)
      real*8, allocatable :: x(:,:),y(:,:),y1(:),gr1(:,:,:),gr2(:,:,:),
     &                       p_mm_0(:,:),p_mm(:,:,:)
      complex*16, allocatable :: xc(:,:),yc(:,:),gc1(:,:,:),gc2(:,:,:)
      ispi1=min(isp1,nspin_1)
      ispi2=min(isp2,nspin_1)
      if(irel<=1) then
        allocate(gr1(nf1,nf2,nd_nrr))
        allocate(gr2(nf1,nf2,nd_nrr))
        call g_mm_sort_r(gr1,1,ind_tau,isp1,nf1,nf2,nd_nrr,n0_nrr,
     &                   nrr_reds,nat1,nat2,isort,jsort,nfun1)
        call g_mm_sort_r(gr2,2,ind_tau,isp2,nf1,nf2,nd_nrr,n0_nrr,
     &                   nrr_reds,nat1,nat2,isort,jsort,nfun1)
      else if(irel==2) then
        allocate(gc1(nf1,nf2,nd_nrr))
        allocate(gc2(nf1,nf2,nd_nrr))
        call g_mm_sort_c(gc1,1,ind_tau,isp1,nf1,nf2,nd_nrr,n0_nrr,
     &                   nrr_reds,nat1,nat2,isort,jsort,nfun1)
        call g_mm_sort_c(gc2,2,ind_tau,isp2,nf1,nf2,nd_nrr,n0_nrr,
     &                   nrr_reds,nat1,nat2,isort,jsort,nfun1)
      endif
      allocate(y1(n2pbm))
      if(irel<=1) then
        allocate(x(nf2,n2pbm),y(nf1,n2pbm))
        allocate(xc(1,1),yc(1,1))
      else if(irel==2) then
        allocate(xc(nf2,n2pbm),yc(nf1,n2pbm))
        allocate(x(1,1),y(1,1))
      endif
      eps=1.d-15
      allocate(nd_p2(nproc_b))
      allocate(n0_p2(nproc_b))
      call size_shift_par(n2pbm,nproc_b,nd_p2,n0_p2)
      allocate(p_mm(n1pbm,nd_p2(me_b+1),nd_nrr))
      p_mm=0.d0
      allocate(p_mm_0(n1pbm,n2pbm))
      ndim0=nf1/nproc_b
      if(ndim0*nproc_b<nf1) ndim0=ndim0+1
      allocate(ndim_km3(nproc_b))
      allocate(index_km3(ndim0,nproc_b))
      do ind_ir0=1,nd_nrr
        ir0=n0_nrr+ind_ir0
        iatom=nrr_reds(2,ir0)
        jatom=nrr_reds(3,ir0)
        ni=nf1
        nipb=n1pbm
        nj=nf2
        njpb=n2pbm
        p_mm_0=0.d0
        call size_shift_par_load(ni,nproc_b,ndim_km3,index_km3,ndim0)
        do ind_km3=1,ndim_km3(me_b+1)
          km3=index_km3(ind_km3,me_b+1)
          lm3=lm_isz(km3,isort)
          lf3=lf_isz(km3,isort)
          if(irel<=1) x=0.d0
          if(irel==2) xc=(0.d0,0.d0)
          do km4=1,nj
            if(irel<=1) g2=gr2(km3,km4,ind_ir0)
            if(irel==2) then
              g2c=conjg(gc2(km3,km4,ind_ir0))
              g2=abs(g2c)
            endif
            if(abs(g2)>=eps) then
              lm4=lm_isz(km4,jsort)
              lf4=lf_isz(km4,jsort)
              do km2=1,nj
                lm2=lm_isz(km2,jsort)
                lf2=lf_isz(km2,jsort)
                if(irel<=1) then
                  ii = max0(lm4,lm2)
                  iii = min0(lm4,lm2)
                  indx = (ii*(ii-1))/2 + iii
                  icg1 = indxcg(indx)
                  icg2 = indxcg(indx+1) - 1
                else
                  icg1=indxcgr(lm2,lm4,1,0)
                  icg2=indxcgr(lm2,lm4,2,0)
                endif
                do icg = icg1, icg2
                  if(irel<=1) l1m1 = jcg(icg)
                  if(irel==2) l1m1 = jcgr(icg)
                  l1=lget(l1m1)
                  if(l1>lmpb(jsort)) cycle
                  if(irel<=1) then
                    cgg=g2*cg(icg)
                    do ie=1,ntle_pb(l1,jsort)
                      lf1=ind_prod0(ie,l1,jsort)
                      km1=ind_prod(ie,l1m1,jsort)
                      x(km2,km1)=x(km2,km1)+cgg
     &                          *fi0(lf2,lf4,lf1,jsort,ispi1,ispi2)
                    enddo
                  else if(irel==2) then
                    cggc=g2c*cgr(icg)
                    do ie=1,ntle_pb(l1,jsort)
                      lf1=ind_prod0(ie,l1,jsort)
                      km1=ind_prod(ie,l1m1,jsort)
                      xc(km2,km1)=xc(km2,km1)+cggc
     &                           *fi0(lf2,lf4,lf1,jsort,ispi1,ispi2)
                    enddo
                  endif
                enddo
              enddo  !! over km2
            endif
          enddo   !! over km4
          if(irel<=1) then
            call dgemm('n','n',nf1,n2pbm,nf2,1.d0,gr1(1,1,ind_ir0),nf1,
     &                 x,nf2,0.d0,y,nf1)
          else if(irel==2) then
            call zgemm('n','n',nf1,n2pbm,nf2,(1.d0,0.d0),
     &                 gc1(1,1,ind_ir0),nf1,xc,nf2,(0.d0,0.d0),yc,nf1)
          endif
          do km1=1,nf1
            lm1=lm_isz(km1,isort)
            lf1=lf_isz(km1,isort)
            if(irel<=1) then
              ii = max0(lm3,lm1)
              iii = min0(lm3,lm1)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1
            else if(irel==2) then
              icg1=indxcgr(lm3,lm1,1,0)
              icg2=indxcgr(lm3,lm1,2,0)
            endif
            do icg = icg1, icg2
              if(irel<=1) lm = jcg(icg)
              if(irel==2) lm = jcgr(icg)
              l=lget(lm)
              if(l>lmpb(isort)) cycle
              do ie=1,ntle_pb(l,isort)
                lf=ind_prod0(ie,l,isort)
                km=ind_prod(ie,lm,isort)
                if(irel<=1) then
                  cgg=fi0(lf1,lf3,lf,isort,ispi1,ispi2)*cg(icg)
                  y1(1:n2pbm)=y(km1,1:n2pbm)*cgg
                else if(irel==2) then
                  cggc=fi0(lf1,lf3,lf,isort,ispi1,ispi2)*cgr(icg)
                  y1(1:n2pbm)=yc(km1,1:n2pbm)*cggc
                endif
                p_mm_0(km,1:n2pbm)=p_mm_0(km,1:n2pbm)-y1(1:n2pbm)
              enddo
            enddo  !! over icg
          enddo  !! over km1
        enddo  !! over ind_km3
        if(nproc_b/=1) call dgop(p_mm_0,n1pbm*n2pbm,'  +',comm_b)
c        ir=nrr_reds(1,ir0)
c       dum=sqrt(dot_product(rr_pnt(:,ir),rr_pnt(:,ir)))*par
c       if(ind_tau==1) print '(2(e12.5,1x))',dum,maxval(abs(p_mm_0))
        do ie=1,nd_p2(me_b+1)
          iii=n0_p2(me_b+1)+ie
          do ii=1,n1pbm
            p_mm(ii,ie,ind_ir0)=p_mm_0(ii,iii)
          enddo
        enddo
      enddo  !! over ind_ir0
      deallocate(ndim_km3,index_km3,y1,p_mm_0,x,y,xc,yc)
      if(irel<=1) then
        deallocate(gr1,gr2)
      else if(irel==2) then
        deallocate(gc1,gc2)
      endif
      call p_mm_sort_fin(p_mm,nd_p2(me_b+1),n0_p2(me_b+1),n1pbm,n2pbm,
     &                   nd_nrr,n0_nrr,isort,jsort,nat1,nat2,nrr_reds,
     &                   r0_pairs,ig_pairs,ind_tau,pmm)
      deallocate(nd_p2,n0_p2,p_mm)
      end


      subroutine p_mm_sclp_0(pwsclp,ind_tau,isp1,isp2,n1pbm,n2pbm,
     &                      nd_nrr,n0_nrr,nfun1,nfun2,isort,jsort,nf1,
     &                      nf2,nrr_maxs,nat1,nat2,nrr_reds,r0_pairs,
     &                      ig_pairs)
      use atom_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      implicit none
      integer, intent(in) :: isp1,isp2,ind_tau,n1pbm,n2pbm,
     &                       nd_nrr,n0_nrr,nfun1,nfun2,isort,jsort,nf1,
     &                       nf2,nrr_maxs,nat1,nat2,
     &                       nrr_reds(3,nqdiv*nat1*nat2),
     &                       r0_pairs(nqdiv,nat1,nat2),
     &                       ig_pairs(nqdiv,nat1,nat2)
      real*8, intent(inout) :: pwsclp(ncmpl*nrow_pw_max,ncol_pw_max,
     &                                ndim3_tn,ndim3_k(me_k+1))
      integer :: lm3,lm,l,ii,iii,indx,lget,icg1,icg2,icg,km3,lf3,
     &           jatom,km2,lm2,lf2,l1m1,l1,km1,lm1,lf1,iatom,
     &           km4,lm4,lf4,ie,lf,km,ir0,ind_ir0,ind_km3,ndim0,ni,
     &           nj,nipb,njpb,ispi1,ispi2
      real*8 :: cgg,eps,g2
      complex*16 :: g2c,cggc
      integer, allocatable :: ndim_km3(:),index_km3(:,:),nd_p2(:),
     &                        n0_p2(:)
      real*8, allocatable :: x(:,:),y(:,:),y1(:),gr1(:,:,:),gr2(:,:,:),
     &                       p_mm_0(:,:),p_mm(:,:,:)
      complex*16, allocatable :: xc(:,:),yc(:,:),gc1(:,:,:),gc2(:,:,:)
      ispi1=min(isp1,nspin_1)
      ispi2=min(isp2,nspin_1)
      if(irel<=1) then
        allocate(gr1(nf1,nf2,nd_nrr))
        allocate(gr2(nf1,nf2,nd_nrr))
        call g_mm_sort_r(gr1,1,ind_tau,isp1,nf1,nf2,nd_nrr,n0_nrr,
     &                   nrr_reds,nat1,nat2,isort,jsort,nfun1)
        call g_mm_sort_r(gr2,2,ind_tau,isp2,nf1,nf2,nd_nrr,n0_nrr,
     &                   nrr_reds,nat1,nat2,isort,jsort,nfun1)
      else if(irel==2) then
        allocate(gc1(nf1,nf2,nd_nrr))
        allocate(gc2(nf1,nf2,nd_nrr))
        call g_mm_sort_c(gc1,1,ind_tau,isp1,nf1,nf2,nd_nrr,n0_nrr,
     &                   nrr_reds,nat1,nat2,isort,jsort,nfun1)
        call g_mm_sort_c(gc2,2,ind_tau,isp2,nf1,nf2,nd_nrr,n0_nrr,
     &                   nrr_reds,nat1,nat2,isort,jsort,nfun1)
      endif
      allocate(y1(n2pbm))
      if(irel<=1) then
        allocate(x(nf2,n2pbm),y(nf1,n2pbm))
        allocate(xc(1,1),yc(1,1))
      else if(irel==2) then
        allocate(xc(nf2,n2pbm),yc(nf1,n2pbm))
        allocate(x(1,1),y(1,1))
      endif
      eps=1.d-15
      allocate(nd_p2(nproc_b))
      allocate(n0_p2(nproc_b))
      call size_shift_par(n2pbm,nproc_b,nd_p2,n0_p2)
      allocate(p_mm(n1pbm,nd_p2(me_b+1),nd_nrr))
      p_mm=0.d0
      allocate(p_mm_0(n1pbm,n2pbm))
      ndim0=nf1/nproc_b
      if(ndim0*nproc_b<nf1) ndim0=ndim0+1
      allocate(ndim_km3(nproc_b))
      allocate(index_km3(ndim0,nproc_b))
      do ind_ir0=1,nd_nrr
        ir0=n0_nrr+ind_ir0
        iatom=nrr_reds(2,ir0)
        jatom=nrr_reds(3,ir0)
        ni=nf1
        nipb=n1pbm
        nj=nf2
        njpb=n2pbm
        p_mm_0=0.d0
        call size_shift_par_load(ni,nproc_b,ndim_km3,index_km3,ndim0)
        do ind_km3=1,ndim_km3(me_b+1)
          km3=index_km3(ind_km3,me_b+1)
          lm3=lm_isz(km3,isort)
          lf3=lf_isz(km3,isort)
          if(irel<=1) x=0.d0
          if(irel==2) xc=(0.d0,0.d0)
          do km4=1,nj
            if(irel<=1) g2=gr2(km3,km4,ind_ir0)
            if(irel==2) then
              g2c=conjg(gc2(km3,km4,ind_ir0))
              g2=abs(g2c)
            endif
            if(abs(g2)>=eps) then
              lm4=lm_isz(km4,jsort)
              lf4=lf_isz(km4,jsort)
              do km2=1,nj
                lm2=lm_isz(km2,jsort)
                lf2=lf_isz(km2,jsort)
                if(irel<=1) then
                  ii = max0(lm4,lm2)
                  iii = min0(lm4,lm2)
                  indx = (ii*(ii-1))/2 + iii
                  icg1 = indxcg(indx)
                  icg2 = indxcg(indx+1) - 1
                else
                  icg1=indxcgr(lm2,lm4,1,0)
                  icg2=indxcgr(lm2,lm4,2,0)
                endif
                do icg = icg1, icg2
                  if(irel<=1) l1m1 = jcg(icg)
                  if(irel==2) l1m1 = jcgr(icg)
                  l1=lget(l1m1)
                  if(l1>lmpb(jsort)) cycle
                  if(irel<=1) then
                    cgg=g2*cg(icg)
                    do ie=1,ntle_pb(l1,jsort)
                      lf1=ind_prod0(ie,l1,jsort)
                      km1=ind_prod(ie,l1m1,jsort)
                      x(km2,km1)=x(km2,km1)+cgg
     &                          *fi0(lf2,lf4,lf1,jsort,ispi1,ispi2)
                    enddo
                  else if(irel==2) then
                    cggc=g2c*cgr(icg)
                    do ie=1,ntle_pb(l1,jsort)
                      lf1=ind_prod0(ie,l1,jsort)
                      km1=ind_prod(ie,l1m1,jsort)
                      xc(km2,km1)=xc(km2,km1)+cggc
     &                           *fi0(lf2,lf4,lf1,jsort,ispi1,ispi2)
                    enddo
                  endif
                enddo
              enddo  !! over km2
            endif
          enddo   !! over km4
          if(irel<=1) then
            call dgemm('n','n',nf1,n2pbm,nf2,1.d0,gr1(1,1,ind_ir0),nf1,
     &                 x,nf2,0.d0,y,nf1)
          else if(irel==2) then
            call zgemm('n','n',nf1,n2pbm,nf2,(1.d0,0.d0),
     &                 gc1(1,1,ind_ir0),nf1,xc,nf2,(0.d0,0.d0),yc,nf1)
          endif
          do km1=1,nf1
            lm1=lm_isz(km1,isort)
            lf1=lf_isz(km1,isort)
            if(irel<=1) then
              ii = max0(lm3,lm1)
              iii = min0(lm3,lm1)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1
            else if(irel==2) then
              icg1=indxcgr(lm3,lm1,1,0)
              icg2=indxcgr(lm3,lm1,2,0)
            endif
            do icg = icg1, icg2
              if(irel<=1) lm = jcg(icg)
              if(irel==2) lm = jcgr(icg)
              l=lget(lm)
              if(l>lmpb(isort)) cycle
              do ie=1,ntle_pb(l,isort)
                lf=ind_prod0(ie,l,isort)
                km=ind_prod(ie,lm,isort)
                if(irel<=1) then
                  cgg=fi0(lf1,lf3,lf,isort,ispi1,ispi2)*cg(icg)
                  y1(1:n2pbm)=y(km1,1:n2pbm)*cgg
                else if(irel==2) then
                  cggc=fi0(lf1,lf3,lf,isort,ispi1,ispi2)*cgr(icg)
                  y1(1:n2pbm)=yc(km1,1:n2pbm)*cggc
                endif
                p_mm_0(km,1:n2pbm)=p_mm_0(km,1:n2pbm)-y1(1:n2pbm)
              enddo
            enddo  !! over icg
          enddo  !! over km1
        enddo  !! over ind_km3
        if(nproc_b/=1) call dgop(p_mm_0,n1pbm*n2pbm,'  +',comm_b)
        do ie=1,nd_p2(me_b+1)
          iii=n0_p2(me_b+1)+ie
          do ii=1,n1pbm
            p_mm(ii,ie,ind_ir0)=p_mm_0(ii,iii)
          enddo
        enddo
      enddo  !! over ind_ir0
      deallocate(ndim_km3,index_km3,y1,p_mm_0,x,y,xc,yc)
      if(irel<=1) then
        deallocate(gr1,gr2)
      else if(irel==2) then
        deallocate(gc1,gc2)
      endif
      call p_mm_sclp_fin(p_mm,nd_p2(me_b+1),n0_p2(me_b+1),n1pbm,n2pbm,
     &                   nd_nrr,n0_nrr,isort,jsort,nat1,nat2,nrr_reds,
     &                   r0_pairs,ig_pairs,ind_tau,pwsclp)
      deallocate(nd_p2,n0_p2,p_mm)
      end
