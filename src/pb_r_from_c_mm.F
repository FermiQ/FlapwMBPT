      subroutine pb_r_from_c_mm(ac,ar,pt)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) :: pt(3)
      complex*16, intent(in) :: ac(n_pbmt,n_pbmt)
      real*8, intent(out) :: ar(n_pbmt,n_pbmt)
      integer :: npw,jcase,icase,jatom,j0,j1,i0,i1,iatom,isort,jsort,l,
     &           l1,jnd0,jnd1,ind0,ind1,ii1,jj1,lm1,lm,nn1,j,i,ii2,jj2,
     &           km,km1,lget
      real*8 :: phase,pi2
      complex*16 :: c1,c2,c3
      complex*16, allocatable :: a(:,:,:)
      pi2=pi+pi
      allocate(a(2,0:maxpb,natom))
      do iatom=1,natom
        isort=is(iatom)
        i0=ip(iatom,inv_num)
        if(i0==iatom) then
          phase=pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))
          do l=0,lmpb(isort)
            a(1,l,iatom)=(0.d0,1.d0)**l*c1
          enddo
        else if(iatom<i0) then
          phase=pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          phase=pi2*dot_product(pt,tau(:,i0))
          c2=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(isort)
            c3=(0.d0,1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*c2
          enddo
        else if(iatom>i0) then
          phase=pi2*dot_product(pt,tau(:,i0))
          c1=-dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          phase=pi2*dot_product(pt,tau(:,iatom))
          c2=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(isort)
            c3=(0.d0,1.d0)**(l+1)
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*c2
          enddo
        endif
      enddo
      do jatom=1,natom
        jsort=is(jatom)
        j0=ip(jatom,inv_num)
        if(j0==jatom) then
          jcase=1
        else if(j0<jatom) then
          jcase=2
          j1=jatom
        else if(j0>jatom) then
          jcase=2
          j1=j0
          j0=jatom
        endif
        jnd0=iopb(j0)-1
        jnd1=jnd0
        if(jcase==2) jnd1=iopb(j1)-1
        do iatom=1,natom
          isort=is(iatom)
          i0=ip(iatom,inv_num)
          if(i0==iatom) then
            icase=1
          else if(i0<iatom) then
            icase=2
            i1=iatom
          else if(i0>iatom) then
            icase=2
            i1=i0
            i0=iatom
          endif
          ind0=iopb(i0)-1
          ind1=ind0
          if(icase==2) ind1=iopb(i1)-1
          do km1=1,n_pbmt0(jsort)
            jj1=jnd0+km1
            jj2=jnd1+km1
            lm1=lm_pbmt(km1,jsort)
            l1=lget(lm1)
            do km=1,n_pbmt0(isort)
              ii1=ind0+km
              ii2=ind1+km
              lm=lm_pbmt(km,isort)
              l=lget(lm)
              if(jcase==1) then
                if(icase==1) then
                  ar(ii1,jj1)=conjg(a(1,l,i0))*a(1,l1,j0)*ac(ii1,jj1)
                else if(icase==2) then
                  ar(ii1,jj1)=(conjg(a(1,l,i0))*ac(ii1,jj1)
     &                        +conjg(a(2,l,i0))*ac(ii2,jj1))*a(1,l1,j0)
                  ar(ii2,jj1)=(conjg(a(1,l,i1))*ac(ii1,jj1)
     &                        +conjg(a(2,l,i1))*ac(ii2,jj1))*a(1,l1,j0)
                endif
              else if(jcase==2) then
                if(icase==1) then
                  ar(ii1,jj1)=conjg(a(1,l,i0))*(a(1,l1,j0)*ac(ii1,jj1)
     &                                         +a(2,l1,j0)*ac(ii1,jj2))
                  ar(ii1,jj2)=conjg(a(1,l,i0))*(a(1,l1,j1)*ac(ii1,jj1)
     &                                         +a(2,l1,j1)*ac(ii1,jj2))
                else if(icase==2) then
                  ar(ii1,jj1)=(conjg(a(1,l,i0))*ac(ii1,jj1)
     &                        +conjg(a(2,l,i0))*ac(ii2,jj1))*a(1,l1,j0)
     &                       +(conjg(a(1,l,i0))*ac(ii1,jj2)
     &                        +conjg(a(2,l,i0))*ac(ii2,jj2))*a(2,l1,j0)
                  ar(ii2,jj1)=(conjg(a(1,l,i1))*ac(ii1,jj1)
     &                        +conjg(a(2,l,i1))*ac(ii2,jj1))*a(1,l1,j0)
     &                       +(conjg(a(1,l,i1))*ac(ii1,jj2)
     &                        +conjg(a(2,l,i1))*ac(ii2,jj2))*a(2,l1,j0)
                  ar(ii1,jj2)=(conjg(a(1,l,i0))*ac(ii1,jj1)
     &                        +conjg(a(2,l,i0))*ac(ii2,jj1))*a(1,l1,j1)
     &                       +(conjg(a(1,l,i0))*ac(ii1,jj2)
     &                        +conjg(a(2,l,i0))*ac(ii2,jj2))*a(2,l1,j1)
                  ar(ii2,jj2)=(conjg(a(1,l,i1))*ac(ii1,jj1)
     &                        +conjg(a(2,l,i1))*ac(ii2,jj1))*a(1,l1,j1)
     &                       +(conjg(a(1,l,i1))*ac(ii1,jj2)
     &                        +conjg(a(2,l,i1))*ac(ii2,jj2))*a(2,l1,j1)
                endif
              endif
            enddo  !! over km
          enddo   !! over km1
        enddo   !! over iatom
      enddo   !! over jatom
      deallocate(a)
      end
