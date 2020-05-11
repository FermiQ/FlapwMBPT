      subroutine p_r_add_mm(pmm,ac,n1pbm,n2pbm,n0,n1,nat1,nat2,isort,
     &                      jsort,pt)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n1pbm,n2pbm,n0,n1,nat1,nat2,isort,jsort
      real*8, intent(in) :: pt(3)
      complex*16, intent(in) :: ac(n1pbm*nat1,n2pbm*nat2)
      real*8, intent(inout) :: pmm(ndim_pbmt(me_b+1))
      integer :: j,km1,m,i,km,ind,ij,iso,i0,i00,i1,i11,j0,j00,j1,j11,
     &           iat,iatom,icase,jat,jatom,jcase,ii1,ii2,ind0,ind1,jj1,
     &           jj2,jnd0,jnd1,l,l1,lm,lm1,lget
      real*8 :: phase,pi2
      complex*16 :: c1,c2,c3
      real*8, allocatable :: ar(:,:)
      complex*16, allocatable :: a(:,:,:)
      pi2=pi+pi
      allocate(a(2,0:maxpb,natom))
      do iatom=1,natom
        iso=is(iatom)
        i0=ip(iatom,inv_num)
        if(i0==iatom) then
          phase=pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))
          do l=0,lmpb(iso)
            a(1,l,iatom)=(0.d0,1.d0)**l*c1
          enddo
        else if(iatom<i0) then
          phase=pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          phase=pi2*dot_product(pt,tau(:,i0))
          c2=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(iso)
            c3=(0.d0,1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*c2
          enddo
        else if(iatom>i0) then
          phase=pi2*dot_product(pt,tau(:,i0))
          c1=-dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          phase=pi2*dot_product(pt,tau(:,iatom))
          c2=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(iso)
            c3=(0.d0,1.d0)**(l+1)
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*c2
          enddo
        endif
      enddo
      allocate(ar(n1pbm*nat1,n2pbm*nat2))
      ar=0.d0
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
      do jat=1,nat2
        jatom=iat_sort(jat,jsort)
        jnd0=iopb(jatom)-1
        j0=n2pbm*(jat-1)
        do iat=1,nat1
          iatom=iat_sort(iat,isort)
          ind0=iopb(iatom)-1
          i0=n1pbm*(iat-1)
          do j=1,n2pbm
            km1=jnd0+j
            m=(km1-1)*km1/2
            do i=1,n1pbm
              km=ind0+i
              if(km>km1) cycle
              ind=m+km
              if(ind<=n0) cycle
              if(ind>n1) cycle
              ij=ind-n0
              pmm(ij)=pmm(ij)+ar(i0+i,j0+j)
            enddo
          enddo
        enddo
      enddo
      deallocate(ar)
      end



      subroutine p_r_add_mm_sclp(pwsclp,ac,n1pbm,n2pbm,n0,n1,nat1,nat2,
     &                           isort,jsort,pt,k)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n1pbm,n2pbm,n0,n1,nat1,nat2,isort,jsort,k
      real*8, intent(in) :: pt(3)
      complex*16, intent(in) :: ac(n1pbm*nat1,n2pbm*nat2)
      real*8, intent(inout) :: pwsclp(nrow_pw_max,ncol_pw_max)
      integer :: j,km1,i,km,iso,i0,i00,i1,i11,j0,j00,j1,j11,
     &           iat,iatom,icase,jat,jatom,jcase,ii1,ii2,ind0,ind1,jj1,
     &           jj2,jnd0,jnd1,l,l1,lm,lm1,lget,ii,jj
      real*8 :: phase,pi2
      complex*16 :: c1,c2,c3
      real*8, allocatable :: ar(:,:)
      complex*16, allocatable :: a(:,:,:)
      pi2=pi+pi
      allocate(a(2,0:maxpb,natom))
      do iatom=1,natom
        iso=is(iatom)
        i0=ip(iatom,inv_num)
        if(i0==iatom) then
          phase=pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))
          do l=0,lmpb(iso)
            a(1,l,iatom)=(0.d0,1.d0)**l*c1
          enddo
        else if(iatom<i0) then
          phase=pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          phase=pi2*dot_product(pt,tau(:,i0))
          c2=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(iso)
            c3=(0.d0,1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*c2
          enddo
        else if(iatom>i0) then
          phase=pi2*dot_product(pt,tau(:,i0))
          c1=-dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          phase=pi2*dot_product(pt,tau(:,iatom))
          c2=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(iso)
            c3=(0.d0,1.d0)**(l+1)
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*c2
          enddo
        endif
      enddo
      allocate(ar(n1pbm*nat1,n2pbm*nat2))
      ar=0.d0
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
      do jat=1,nat2
        jatom=iat_sort(jat,jsort)
        jnd0=iopb(jatom)-1
        j0=n2pbm*(jat-1)
        do iat=1,nat1
          iatom=iat_sort(iat,isort)
          ind0=iopb(iatom)-1
          i0=n1pbm*(iat-1)
          do j=1,ncol_pw(k)
            jj=refcol_pw(j,k)-jnd0
            if(jj<1.or.jj>n2pbm) cycle
            do i=1,nrow_pw(k)
              ii=refrow_pw(i,k)-ind0
              if(ii<1.or.ii>n1pbm) cycle
              pwsclp(i,j)=pwsclp(i,j)+ar(i0+ii,j0+jj)
            enddo
          enddo
          if(isort/=jsort) then
            do j=1,ncol_pw(k)
              jj=refcol_pw(j,k)-ind0
              if(jj<1.or.jj>n1pbm) cycle
              do i=1,nrow_pw(k)
                ii=refrow_pw(i,k)-jnd0
                if(ii<1.or.ii>n2pbm) cycle
                pwsclp(i,j)=pwsclp(i,j)+ar(i0+jj,j0+ii)
              enddo
            enddo
          endif
        enddo
      enddo
      deallocate(ar)
      end
