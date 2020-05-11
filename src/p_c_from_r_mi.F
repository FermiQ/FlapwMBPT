      subroutine p_c_from_r_mi(npw,ac,ar,k)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: npw,k
      real*8, intent(in) :: ar(nd_b_pbmt(me_b+1),npw)
      complex*16, intent(out) :: ac(n_pbmt,npw)
      integer :: icase,i0,i1,iatom,isort,l,ind0,ind1,ii1,lm,ii2,km,lget,
     &           n0,n1,i
      real*8 :: phase,pi2,pt(3)
      complex*16 :: c1,c3
      complex*16, allocatable :: a(:,:,:)
      pi2=pi+pi
      n0=n0_b_pbmt(me_b+1)
      n1=n0+nd_b_pbmt(me_b+1)
      pt=pnt(:,k)
      allocate(a(2,0:maxpb,natom))
      do iatom=1,natom
        isort=is(iatom)
        i0=ip(iatom,inv_num)
        if(i0==iatom) then
          phase=-pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))
          do l=0,lmpb(isort)
            a(1,l,iatom)=(0.d0,-1.d0)**l*c1
          enddo
        else if(i0<iatom) then
          phase=-pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(isort)
            c3=(0.d0,-1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*(0.d0,-1.d0)*c1
          enddo
        else if(i0>iatom) then
          phase=-pi2*dot_product(pt,tau(:,iatom))
          c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
          do l=0,lmpb(isort)
            c3=(0.d0,-1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=-c3*(0.d0,-1.d0)*c1
          enddo
        endif
      enddo
      ac=(0.d0,0.d0)
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
            if(ii1>n0.and.ii1<=n1) then
              i=ii1-n0
              ac(ii1,:)=conjg(a(1,l,i0))*ar(i,:)
            endif
          else if(icase==2) then
            if(ii1>n0.and.ii1<=n1) then
              i=ii1-n0
              ac(ii1,:)=ac(ii1,:)+conjg(a(1,l,i0))*ar(i,:)
              ac(ii2,:)=ac(ii2,:)+conjg(a(1,l,i1))*ar(i,:)
            endif
            if(ii2>n0.and.ii2<=n1) then
              i=ii2-n0
              ac(ii1,:)=ac(ii1,:)+conjg(a(2,l,i0))*ar(i,:)
              ac(ii2,:)=ac(ii2,:)+conjg(a(2,l,i1))*ar(i,:)
            endif
          endif
        enddo   !! over km
      enddo   !! over iatom
      deallocate(a)
      if(nproc_b/=1) call dgop(ac,2*n_pbmt*npw,'  +',comm_b)
      end
