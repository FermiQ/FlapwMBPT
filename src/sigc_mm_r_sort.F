      subroutine sigc_mm_r_sort(ispin,it,ind_tau,s_c,w_rs,n1pbm,n2pbm,
     &                          nd_nrr,n0_nrr,nfun1,nfun2,isort,jsort,
     &                          nf1,nf2,nrr_maxs,nat1,nat2,nrr_reds,
     &                          r0_pairs,ig_pairs)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ispin,n1pbm,n2pbm,nd_nrr,n0_nrr,nfun1,
     &                       nfun2,isort,jsort,nf1,nf2,nrr_maxs,nat1,
     &                       nat2,nrr_reds(3,nqdiv*nat1*nat2),it,
     &                       r0_pairs(nqdiv,nat1,nat2),ind_tau,
     &                       ig_pairs(nqdiv,nat1,nat2)
      real*8, intent(in) :: w_rs(n1pbm,n2pbm,nd_nrr)
      complex*16, intent(out) :: s_c(nfun1,nfun2,ndim3_kk(me_k+1))
      integer :: lm3,lm,ii,iii,indx,lget,icg1,icg2,icg,km3,lf3,jatom,
     &           km2,lm2,lf2,km1,lm1,lf1,iatom,ir0,ir,lm22,l22,km,lf,
     &           lm11,l11,ind_ir0,isp,ie,lf22,km22,lf11,km11,ind_km1,
     &           ndim0
      real*8 :: cgg,dum
      complex*16 :: cggc,dumc
      integer, allocatable :: ndim_km1(:),index_km1(:,:)
      real*8, allocatable :: x(:,:),y(:,:),sr(:,:,:),gr(:,:,:)
      complex*16, allocatable :: xc(:,:),yc(:,:),sc(:,:,:),vtmp(:,:),
     &                         gc(:,:,:)
      isp=min(ispin,nspin_1)
c ---- Prepare Green's functions in real space ---------------------
      if(irel<=1) then
        allocate(gr(nf1,nf2,nd_nrr))
        call g_rs_from_ks_sort(gr,dumc,it,ind_tau,ispin,nf1,nf2,nd_nrr,
     &                         n0_nrr,nfun1,nfun2,isort,jsort,nrr_maxs,
     &                         nat1,nat2,nrr_reds)
      else if(irel==2) then
        allocate(gc(nf1,nf2,nd_nrr))
        call g_rs_from_ks_sort(dum,gc,it,ind_tau,ispin,nf1,nf2,nd_nrr,
     &                         n0_nrr,nfun1,nfun2,isort,jsort,nrr_maxs,
     &                         nat1,nat2,nrr_reds)
      endif
      if(irel<=1) then
        allocate(sr(nf1,nf2,nrr_maxs))
        allocate(x(n2pbm,nf1),y(n1pbm,nf1))
        sr=0.d0
      else if(irel==2) then
        allocate(sc(nf1,nf2,nrr_maxs))
        allocate(xc(n2pbm,nf1),yc(n1pbm,nf1))
        sc=(0.d0,0.d0)
        allocate(vtmp(n1pbm,n2pbm))
      endif
      ndim0=nf2/nproc_b
      if(ndim0*nproc_b<nf2) ndim0=ndim0+1
      allocate(ndim_km1(nproc_b))
      allocate(index_km1(ndim0,nproc_b))
      do ind_ir0=1,nd_nrr
        ir0=n0_nrr+ind_ir0
        ir=nrr_reds(1,ir0)
        iatom=nrr_reds(2,ir0)
        jatom=nrr_reds(3,ir0)
        if(irel==2) vtmp=w_rs(:,:,ind_ir0)
        call size_shift_par_load(nf2,nproc_b,ndim_km1,index_km1,ndim0)
        do ind_km1=1,ndim_km1(me_b+1)
          km1=index_km1(ind_km1,me_b+1)
          lm1=lm_isz(km1,jsort)
          lf1=lf_isz(km1,jsort)
          if(irel<=1) x=0.d0
          if(irel==2) xc=(0.d0,0.d0)
          do km3=1,nf2
            lm3=lm_isz(km3,jsort)
            lf3=lf_isz(km3,jsort)
            if(irel<=1) then
              ii = max0(lm3,lm1)
              iii = min0(lm3,lm1)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1
            else
              icg1=indxcgr(lm3,lm1,1,0)
              icg2=indxcgr(lm3,lm1,2,0)
            endif
            do icg = icg1, icg2
              if(irel<=1) lm22 = jcg(icg)
              if(irel==2) lm22 = jcgr(icg)
              l22=lget(lm22)
              if(l22>lmpb(jsort)) cycle
              do ie=1,ntle_pb(l22,jsort)
                lf22=ind_prod0(ie,l22,jsort)
                km22=ind_prod(ie,lm22,jsort)
                if(irel<=1) then
                  cgg=fi0(lf3,lf1,lf22,jsort,isp,isp)*cg(icg)
                  x(km22,:)=x(km22,:)+gr(:,km3,ind_ir0)*cgg
                else if(irel==2) then
                  cggc=fi0(lf3,lf1,lf22,jsort,isp,isp)*cgr(icg)
                  xc(km22,:)=xc(km22,:)+gc(:,km3,ind_ir0)*cggc
                endif
              enddo   !! over ie
            enddo   !! over icg
          enddo   !! over km3
          if(irel<=1) then
            call dgemm('n','n',n1pbm,nf1,n2pbm,1.d0,w_rs(1,1,ind_ir0),
     &                 n1pbm,x,n2pbm,0.d0,y,n1pbm)
          else if(irel==2) then
            call zgemm('n','n',n1pbm,nf1,n2pbm,(1.d0,0.d0),vtmp,n1pbm,
     &                 xc,n2pbm,(0.d0,0.d0),yc,n1pbm)
          endif
          do km=1,nf1
            lm=lm_isz(km,isort)
            lf=lf_isz(km,isort)
            do km2=1,nf1
              lm2=lm_isz(km2,isort)
              lf2=lf_isz(km2,isort)
              if(irel<=1) then
                ii = max0(lm2,lm)
                iii = min0(lm2,lm)
                indx = (ii*(ii-1))/2 + iii
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1) - 1
              else
                icg1=indxcgr(lm,lm2,1,0)
                icg2=indxcgr(lm,lm2,2,0)
              endif
              do icg = icg1, icg2
                if(irel<=1) lm11 = jcg(icg)
                if(irel==2) lm11 = jcgr(icg)
                l11=lget(lm11)
                if(l11>lmpb(isort)) cycle
                do ie=1,ntle_pb(l11,isort)
                  lf11=ind_prod0(ie,l11,isort)
                  km11=ind_prod(ie,lm11,isort)
                  if(irel<=1) then
                    cgg=fi0(lf2,lf,lf11,isort,isp,isp)*cg(icg)
                    sr(km,km1,ir0)=sr(km,km1,ir0)-cgg*y(km11,km2)
                  else if(irel==2) then
                    cggc=fi0(lf2,lf,lf11,isort,isp,isp)*cgr(icg)
                    sc(km,km1,ir0)=sc(km,km1,ir0)-cggc*yc(km11,km2)
                  endif
                enddo
              enddo
            enddo
         enddo  !! over km
       enddo    !! over ind_km1
c       dum=sqrt(dot_product(rr_pnt(:,ir),rr_pnt(:,ir)))*par
c       if(it==1.and.ind_tau==1) print '(2(e12.5,1x))',dum,
c     &                           maxval(abs(sr(:,:,ir0)))
      enddo  !! over ind_ir0
      deallocate(ndim_km1,index_km1)
      if(nproc_b/=1) then
        if(irel<=1) call dgop(sr,nf1*nf2*nrr_maxs,'  +',comm_b)
        if(irel==2) call dgop(sc,2*nf1*nf2*nrr_maxs,'  +',comm_b)
      endif
      if(nproc_k/=1) then
        if(irel<=1) call dgop(sr,nf1*nf2*nrr_maxs,'  +',comm_k)
        if(irel==2) call dgop(sc,2*nf1*nf2*nrr_maxs,'  +',comm_k)
      endif
      if(irel<=1) then
        deallocate(x,y,gr)
      else if(irel==2) then
        deallocate(xc,yc,vtmp,gc)
      endif
      if(irel<=1) then
        call sym_sig_rs_sort(sr,s_c,nf1,nf2,nfun1,nfun2,isort,jsort,
     &                       nrr_maxs,nat1,nat2,r0_pairs,ig_pairs)
        deallocate(sr)
      else if(irel==2) then
        call sym_sig_rs_rel_sort(sc,s_c,nf1,nf2,nfun1,nfun2,isort,jsort,
     &                           nrr_maxs,nat1,nat2,r0_pairs,ig_pairs)
        deallocate(sc)
      endif
      end
