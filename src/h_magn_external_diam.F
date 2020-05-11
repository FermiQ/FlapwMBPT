      subroutine h_magn_external_diam
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      IMPLICIT none
      INTEGER :: lf1,mt1,lf,mt,irad,isort,km1,lm1,l1,i1,mj1,li1,jj1,km,
     &           lm,l,i,mj,li,jj,iatom,lmm,lmm1,lm2,icg,icg1,icg2
      real*8 :: work(0:maxnrad),f11,f22,c2,dqdall,b(3),d1,d2,d3,d4,const
      complex*16 :: cc,xx,yy,zz,xy,xz,yz
      real*8, allocatable :: f1(:,:),f2(:,:)
      allocate(f1(maxlfun,maxlfun),f2(maxlfun,maxlfun))
      c2=clight**2
      const=2.d0*c2
      d1=sqrt(16.d0*pi/9.d0)
      d2=sqrt(4.d0*pi/15.d0)
      d3=sqrt(4.d0*pi/45.d0)
      d4=sqrt(16.d0*pi/45.d0)
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
          call getlimj(lmm1,l1+i1,-i1,mj1,li1,1)
          jj1=l1+l1+i1
          do km=1,lfunm(isort)
            lm=lm_isz(km,isort)
            lf=lf_isz(km,isort)
            call getlimj(lm,l,i,mj,li,0)
            call getlimj(lmm,l+i,-i,mj,li,1)
            jj=l+l+i
            f11=f1(lf,lf1)
            f22=f2(lf,lf1)
            xx=(0.d0,0.d0)
            yy=(0.d0,0.d0)
            zz=(0.d0,0.d0)
            xy=(0.d0,0.d0)
            xz=(0.d0,0.d0)
            yz=(0.d0,0.d0)
c -------- Big components -------------------------
            icg1=indxcgr(lm,lm1,1,0)
            icg2=indxcgr(lm,lm1,2,0)
            do icg=icg1,icg2
              lm2=jcgr(icg)
              cc=cgr(icg)
              if(lm2==1) then
                xx=xx+d1*cc*f11
                yy=yy+d1*cc*f11
                zz=zz+d1*cc*f11
              else if(lm2==5) then
                xy=xy+d2*cc*f11
              else if(lm2==6) then
                yz=yz+d2*cc*f11
              else if(lm2==7) then
                xx=xx+d3*cc*f11
                yy=yy+d3*cc*f11
                zz=zz-d4*cc*f11
              else if(lm2==8) then
                xz=xz+d2*cc*f11
              else if(lm2==9) then
                xx=xx-d2*cc*f11
                yy=yy+d2*cc*f11
              endif
            enddo
c -------- Small components -------------------------
            icg1=indxcgr(lmm,lmm1,1,0)
            icg2=indxcgr(lmm,lmm1,2,0)
            do icg=icg1,icg2
              lm2=jcgr(icg)
              cc=cgr(icg)
              if(lm2==1) then
                xx=xx+d1*cc*f22
                yy=yy+d1*cc*f22
                zz=zz+d1*cc*f22
              else if(lm2==5) then
                xy=xy+d2*cc*f22
              else if(lm2==6) then
                yz=yz+d2*cc*f22
              else if(lm2==7) then
                xx=xx+d3*cc*f22
                yy=yy+d3*cc*f22
                zz=zz-d4*cc*f22
              else if(lm2==8) then
                xz=xz+d2*cc*f22
              else if(lm2==9) then
                xx=xx-d2*cc*f22
                yy=yy+d2*cc*f22
              endif
            enddo
            cc=xx*b(1)**2+yy*b(2)**2+zz*b(3)**2-2.d0*xy*b(1)*b(2)
     &        -2.d0*xz*b(1)*b(3)-2.d0*yz*b(2)*b(3)
            pv(km,km1,iatom,1)=pv(km,km1,iatom,1)+real(cc)/const
            pvj(km,km1,iatom)=pvj(km,km1,iatom)+imag(cc)/const
          enddo  !!! over km
        enddo  !!! over km1
      enddo   !! over iatom
      deallocate(f1,f2)
      end
