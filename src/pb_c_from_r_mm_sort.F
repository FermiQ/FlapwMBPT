      subroutine pb_c_from_r_mm_sort(ac,ar,pt,isort,jsort,nat1,nat2,
     &                               n1pbm,n2pbm)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: isort,jsort,nat1,nat2,n1pbm,n2pbm
      real*8, intent(in) :: pt(3),ar(n1pbm*nat1,n2pbm*nat2)
      complex*16, intent(out) :: ac(n1pbm*nat1,n2pbm*nat2)
      integer :: jcase,icase,jatom,j0,j1,i0,i1,iatom,l,l1,i00,i11,j00,
     &           j11,jnd0,jnd1,ind0,ind1,ii1,jj1,lm1,lm,ii2,jj2,km,km1,
     &           lget,iat,jat,iso
      real*8 :: phase,pi2
      complex*16 :: c1,c3
      complex*16, allocatable :: a(:,:,:)
      pi2=pi+pi
      allocate(a(2,0:maxpb,natom))
      do iatom=1,natom
        iso=is(iatom)
        i0=ip(iatom,inv_num)
        if(i0==iatom) then
          phase=-pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))
          do l=0,lmpb(iso)
            a(1,l,iatom)=(0.d0,-1.d0)**l*c1
          enddo
        else if(i0<iatom) then
          phase=-pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(iso)
            c3=(0.d0,-1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*(0.d0,-1.d0)*c1
          enddo
        else if(i0>iatom) then
          phase=-pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(iso)
            c3=(0.d0,-1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=-c3*(0.d0,-1.d0)*c1
          enddo
        endif
      enddo
      do jat=1,nat2
        jatom=iat_sort(jat,jsort)
        j0=ip(jatom,inv_num)
        if(j0==jatom) then
          jcase=1
          j1=j0
        else if(j0<jatom) then
          jcase=2
          j1=jatom
        else if(j0>jatom) then
          jcase=2
          j1=j0
          j0=jatom
        endif
        j00=iat_sort_back(j0)
        j11=iat_sort_back(j1)
        jnd0=n2pbm*(j00-1)
        jnd1=jnd0
        if(jcase==2) jnd1=n2pbm*(j11-1)
        do iat=1,nat1
          iatom=iat_sort(iat,isort)
          i0=ip(iatom,inv_num)
          if(i0==iatom) then
            icase=1
            i1=i0
          else if(i0<iatom) then
            icase=2
            i1=iatom
          else if(i0>iatom) then
            icase=2
            i1=i0
            i0=iatom
          endif
          i00=iat_sort_back(i0)
          i11=iat_sort_back(i1)
          ind0=n1pbm*(i00-1)
          ind1=ind0
          if(icase==2) ind1=n1pbm*(i11-1)
          do km1=1,n2pbm
            jj1=jnd0+km1
            jj2=jnd1+km1
            lm1=lm_pbmt(km1,jsort)
            l1=lget(lm1)
            do km=1,n1pbm
              ii1=ind0+km
              ii2=ind1+km
              lm=lm_pbmt(km,isort)
              l=lget(lm)
              if(jcase==1) then
                if(icase==1) then
                  ac(ii1,jj1)=conjg(a(1,l,i0))*a(1,l1,j0)*ar(ii1,jj1)
                else if(icase==2) then
                  ac(ii1,jj1)=(conjg(a(1,l,i0))*ar(ii1,jj1)
     &                        +conjg(a(2,l,i0))*ar(ii2,jj1))*a(1,l1,j0)
                  ac(ii2,jj1)=(conjg(a(1,l,i1))*ar(ii1,jj1)
     &                        +conjg(a(2,l,i1))*ar(ii2,jj1))*a(1,l1,j0)
                endif
              else if(jcase==2) then
                if(icase==1) then
                  ac(ii1,jj1)=conjg(a(1,l,i0))*(a(1,l1,j0)*ar(ii1,jj1)
     &                                         +a(2,l1,j0)*ar(ii1,jj2))
                  ac(ii1,jj2)=conjg(a(1,l,i0))*(a(1,l1,j1)*ar(ii1,jj1)
     &                                         +a(2,l1,j1)*ar(ii1,jj2))
                else if(icase==2) then
                  ac(ii1,jj1)=(conjg(a(1,l,i0))*ar(ii1,jj1)
     &                        +conjg(a(2,l,i0))*ar(ii2,jj1))*a(1,l1,j0)
     &                       +(conjg(a(1,l,i0))*ar(ii1,jj2)
     &                        +conjg(a(2,l,i0))*ar(ii2,jj2))*a(2,l1,j0)
                  ac(ii2,jj1)=(conjg(a(1,l,i1))*ar(ii1,jj1)
     &                        +conjg(a(2,l,i1))*ar(ii2,jj1))*a(1,l1,j0)
     &                       +(conjg(a(1,l,i1))*ar(ii1,jj2)
     &                        +conjg(a(2,l,i1))*ar(ii2,jj2))*a(2,l1,j0)
                  ac(ii1,jj2)=(conjg(a(1,l,i0))*ar(ii1,jj1)
     &                        +conjg(a(2,l,i0))*ar(ii2,jj1))*a(1,l1,j1)
     &                       +(conjg(a(1,l,i0))*ar(ii1,jj2)
     &                        +conjg(a(2,l,i0))*ar(ii2,jj2))*a(2,l1,j1)
                  ac(ii2,jj2)=(conjg(a(1,l,i1))*ar(ii1,jj1)
     &                        +conjg(a(2,l,i1))*ar(ii2,jj1))*a(1,l1,j1)
     &                       +(conjg(a(1,l,i1))*ar(ii1,jj2)
     &                        +conjg(a(2,l,i1))*ar(ii2,jj2))*a(2,l1,j1)
                endif
              endif
            enddo  !! over km
          enddo   !! over km1
        enddo   !! over iat
      enddo   !! over jat
      deallocate(a)
      end
