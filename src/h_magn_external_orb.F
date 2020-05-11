      subroutine h_magn_external_orb
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      IMPLICIT none
      INTEGER :: lf1,mt1,lf,mt,irad,isort,km1,lm1,l1,i1,mj1,li1,jj1,km,
     &           lm,l,i,mj,li,jj,iatom
      real*8 :: work(0:maxnrad),f11,f22,c2,dqdall,aorbx,aorby,aorbz,
     &          aorbxi,aorbyi,aorbzi,b(3),fc
      real*8, allocatable :: f1(:,:),f2(:,:)
      allocate(f1(maxlfun,maxlfun),f2(maxlfun,maxlfun))
      c2=clight**2
      b=b_ext*b_extval
      do iatom=1,natom
        isort=is(iatom)
        do lf1=1,lfun(isort)
          mt1=ind_wf(lf1,isort)
          do lf=1,lfun(isort)
            mt=ind_wf(lf,isort)
            do irad=0,nrad(isort)
              work(irad)=gfun(mt+irad,1)*gfun(mt1+irad,1)*dr(irad,isort)
     &                                           *r(irad,isort)**2
            enddo
            f1(lf,lf1)=dqdall(h(isort),work,nrad(isort))
            do irad=0,nrad(isort)
              work(irad)=-gfund(mt+irad,1)*gfund(mt1+irad,1)
     &                   *dr(irad,isort)*r(irad,isort)**2
            enddo
            f2(lf,lf1)=dqdall(h(isort),work,nrad(isort))/c2
          enddo  !!! over lf
        enddo  !!! over lf1
        do km1=1,lfunm(isort)
          lm1=lm_isz(km1,isort)
          lf1=lf_isz(km1,isort)
          call getlimj(lm1,l1,i1,mj1,li1,0)
          jj1=l1+l1+i1
          do km=1,lfunm(isort)
            lm=lm_isz(km,isort)
            lf=lf_isz(km,isort)
            call getlimj(lm,l,i,mj,li,0)
            jj=l+l+i
            f11=f1(lf,lf1)
            f22=f2(lf,lf1)
            call orb_rel(l,i,mj,l1,i1,mj1,aorbx,aorby,aorbz,aorbxi,
     &                   aorbyi,aorbzi)
            fc=f11*(b(1)*aorbx+b(2)*aorby+b(3)*aorbz)
            pv(km,km1,iatom,1)=pv(km,km1,iatom,1)+fc
            fc=f11*(b(1)*aorbxi+b(2)*aorbyi+b(3)*aorbzi)
            pvj(km,km1,iatom)=pvj(km,km1,iatom)+fc
            call orb_rel(l+i,-i,mj,l1+i1,-i1,mj1,aorbx,aorby,aorbz,
     &                   aorbxi,aorbyi,aorbzi)
            fc=f22*(b(1)*aorbx+b(2)*aorby+b(3)*aorbz)
            pv(km,km1,iatom,1)=pv(km,km1,iatom,1)+fc
            fc=f22*(b(1)*aorbxi+b(2)*aorbyi+b(3)*aorbzi)
            pvj(km,km1,iatom)=pvj(km,km1,iatom)+fc
          enddo  !!! over km
        enddo  !!! over km1
      enddo   !! over iatom
      deallocate(f1,f2)
      end
