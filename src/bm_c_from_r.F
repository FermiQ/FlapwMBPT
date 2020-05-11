      subroutine bm_c_from_r(npw,ac,ar)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: npw
      real*8, intent(in) :: ar(n_opt_pb,n_pbmt+npw)
      complex*16, intent(out) :: ac(n_opt_pb,n_pbmt+npw)
      integer :: icase,i0,i1,iatom,isort,l,ind0,ind1,ii1,lm,ii2,km,lget,
     &           j
      complex*16 :: c1,c3
      complex*16, allocatable :: a(:,:,:)
      allocate(a(2,0:maxpb,natom))
      do iatom=1,natom
        isort=is(iatom)
        i0=ip(iatom,inv_num)
        if(i0==iatom) then
          do l=0,lmpb(isort)
            a(1,l,iatom)=(0.d0,-1.d0)**l
          enddo
        else if(i0<iatom) then
          c1=1.d0/sqrt(2.d0)
          do l=0,lmpb(isort)
            c3=(0.d0,-1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=c3*(0.d0,-1.d0)*c1
          enddo
        else if(i0>iatom) then
          c1=1.d0/sqrt(2.d0)
          do l=0,lmpb(isort)
            c3=(0.d0,-1.d0)**l
            a(1,l,iatom)=c3*c1
            a(2,l,iatom)=-c3*(0.d0,-1.d0)*c1
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
            ac(:,ii1)=a(1,l,i0)*ar(:,ii1)
          else if(icase==2) then
            ac(:,ii1)=a(1,l,i0)*ar(:,ii1)+a(2,l,i0)*ar(:,ii2)
            ac(:,ii2)=a(1,l,i1)*ar(:,ii1)+a(2,l,i1)*ar(:,ii2)
          endif
        enddo   !! over km
      enddo   !! over iatom
      deallocate(a)
      do j=n_pbmt+1,n_pbmt+npw
        ac(:,j)=ar(:,j)
      enddo
      end
