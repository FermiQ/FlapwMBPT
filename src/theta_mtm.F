      subroutine theta_mtm
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: isort,iatom,nlmf,j0,jfn,i0,ifn,lm2,isym,mt,ii,
     &	       irad,j,lmj,i,lmi,ix,icg1,icg2,icg
      real*8 :: cgg,dqdall,sqpi4,v0
      real*8, allocatable :: t(:,:,:,:),work(:),work0(:)
      thet_me=0.d0
      allocate(t(maxlfpb,maxlfpb,maxnsym,n_ixc_0),work(0:maxnrad))
      allocate(work0(0:maxnrad))
      sqpi4=sqrt(4.d0*pi)
      do isort=1,nsort
        iatom=iat_1(isort)
        nlmf=(lmpb(isort)+1)**2
        do j0=1,lfun_pb(isort)
          jfn=ind_pb(j0,isort)
          do i0=1,lfun_pb(isort)
            ifn=ind_pb(i0,isort)
            do lm2=1,nlmf
              if(sym(lm2,isort)) then
                isym=lmsym(lm2,isort)
                mt=indmt(isym,isort,1)
                do irad=0,nrad(isort)
                  work0(irad)=fun_pb(irad+ifn)*fun_pb(irad+jfn)
                enddo
                do i=1,n_ixc_0
                  do irad=0,nrad(isort)
                    v0=theta_mt(mt+irad,i)
                    work(irad)=work0(irad)*v0*dr(irad,isort)
     &			                             *r(irad,isort)**2
                  enddo
                  t(i0,j0,isym,i)=dqdall(h(isort),work,nrad(isort))
                enddo
              endif
            enddo   !! over lm2
          enddo   !! over i0
        enddo   !! over j0
        do j=1,n_pbmt0(isort)
          j0=lf_pb(j,isort)
          lmj=lm_pbmt(j,isort)
          do i=1,n_pbmt0(isort)
            i0=lf_pb(i,isort)
            lmi=lm_pbmt(i,isort)
            ix=max0(lmi,lmj)
            ix=(ix*(ix-1))/2+min0(lmi,lmj)
            icg1=indxcg(ix)
            icg2=indxcg(ix+1)-1
            do icg=icg1,icg2
              lm2 = jcg(icg)
              cgg=cg(icg)
              if(lm2.le.nlmf) then
                if(sym(lm2,isort)) then
                  isym=lmsym(lm2,isort)
                  do ii=1,n_ixc_0
                    thet_me(i,j,iatom,ii)=thet_me(i,j,iatom,ii)
     &			                       +cgg*t(i0,j0,isym,ii)
                  enddo
                endif
              endif
            enddo  !! over icg
          enddo   !! over i
        enddo   !! over j
      enddo  !! over isort
      deallocate(t,work,work0)
c ----------------------------------------------------------------	  
      call sym_me_pb
      end