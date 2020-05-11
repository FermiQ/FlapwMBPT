      subroutine p_r_from_c_mi(npw,ac,ar,k)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: npw,k
      complex*16, intent(in) :: ac(nr_full,n_pbmt)
      real*8, intent(out) :: ar(n_pbmt,npw)
      integer :: icase,i0,i1,iatom,isort,l,ind0,ind1,ii1,lm,j,ii2,
     &           km,lget,ind,j0
      real*8 :: phase,pi2,pt(3)
      complex*16 :: c1,c2,c3
      complex*16, allocatable :: a(:,:,:)
      pi2=pi+pi
      pt=pnt(:,k)
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
      do j=1,npw
        j0=indpw_gw(j,k)
        call zone1_number(gbs(:,j0),rbas,nrdiv,ind)
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
          do km=1,n_pbmt0(isort)
            ii1=ind0+km
            ii2=ind1+km
            lm=lm_pbmt(km,isort)
            l=lget(lm)
            if(icase==1) then
              ar(ii1,j)=conjg(a(1,l,i0))*ac(ind,ii1)
            else if(icase==2) then
              ar(ii1,j)=conjg(a(1,l,i0))*ac(ind,ii1)
     &                 +conjg(a(2,l,i0))*ac(ind,ii2)
              ar(ii2,j)=conjg(a(1,l,i1))*ac(ind,ii1)
     &                 +conjg(a(2,l,i1))*ac(ind,ii2)
            endif
          enddo   !! over km
        enddo   !! over iatom
      enddo
      deallocate(a)
      end


      subroutine p_r_from_c_mi_0(npw,ac,ar,k)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: npw,k
      complex*16, intent(in) :: ac(n_pbmt,npw)
      real*8, intent(out) :: ar(n_pbmt,npw)
      integer :: icase,i0,i1,iatom,isort,l,ind0,ind1,ii1,lm,ii2,km,lget
      real*8 :: phase,pi2,pt(3)
      complex*16 :: c1,c2,c3
      complex*16, allocatable :: a(:,:,:)
      pi2=pi+pi
      pt=pnt(:,k)
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
        do km=1,n_pbmt0(isort)
          ii1=ind0+km
          ii2=ind1+km
          lm=lm_pbmt(km,isort)
          l=lget(lm)
          if(icase==1) then
            ar(ii1,:)=conjg(a(1,l,i0))*ac(ii1,:)
          else if(icase==2) then
            ar(ii1,:)=conjg(a(1,l,i0))*ac(ii1,:)
     &               +conjg(a(2,l,i0))*ac(ii2,:)
            ar(ii2,:)=conjg(a(1,l,i1))*ac(ii1,:)
     &               +conjg(a(2,l,i1))*ac(ii2,:)
          endif
        enddo   !! over km
      enddo   !! over iatom
      deallocate(a)
      end
