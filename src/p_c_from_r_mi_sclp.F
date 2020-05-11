      subroutine p_c_from_r_mi_sclp(npw,ac,ar,k)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: npw,k
      real*8, intent(in) :: ar(nrow_pw_max,ncol_pw_max)
      complex*16, intent(out) :: ac(n_pbmt,npw)
      integer :: icase,i0,i1,iatom,isort,l,ind0,ind1,ii1,lm,ii2,km,lget,
     &           i,j,jj
      real*8 :: phase,pi2,pt(3)
      complex*16 :: c1,c3
      integer, allocatable :: km_row(:)
      complex*16, allocatable :: a(:,:,:)
      pi2=pi+pi
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
      allocate(km_row(n_pbmt))
      km_row=0
      do i=1,nrow_pw(k)
        j=refrow_pw(i,k)
        if(j<=n_pbmt) km_row(j)=i
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
            i=km_row(ii1)
            if(i/=0) then
              do j=1,ncol_pw(k)
                jj=refcol_pw(j,k)-n_pbmt
                if(jj>0) ac(ii1,jj)=conjg(a(1,l,i0))*ar(i,j)
              enddo
            endif
          else if(icase==2) then
            i=km_row(ii1)
            if(i/=0) then
              do j=1,ncol_pw(k)
                jj=refcol_pw(j,k)-n_pbmt
                if(jj<=0) cycle
                ac(ii1,jj)=ac(ii1,jj)+conjg(a(1,l,i0))*ar(i,j)
                ac(ii2,jj)=ac(ii2,jj)+conjg(a(1,l,i1))*ar(i,j)
              enddo
            endif
            i=km_row(ii2)
            if(i/=0) then
              do j=1,ncol_pw(k)
                jj=refcol_pw(j,k)-n_pbmt
                if(jj<=0) cycle
                ac(ii1,jj)=ac(ii1,jj)+conjg(a(2,l,i0))*ar(i,j)
                ac(ii2,jj)=ac(ii2,jj)+conjg(a(2,l,i1))*ar(i,j)
              enddo
            endif
          endif
        enddo   !! over km
      enddo   !! over iatom
      deallocate(a,km_row)
      call dgop(ac,2*n_pbmt*npw,'  +',comm_b)
      end
