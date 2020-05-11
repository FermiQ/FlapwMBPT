      subroutine theta_r_from_c(ac,ar)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) :: ac(n_pbmtm,n_pbmtm,natom)
      real*8, intent(out) :: ar(n_pbmtm,n_pbmtm,natom)
      integer :: icase,i0,i1,iatom,isort,l,l1,ind0,ind1,lm1,lm,km,km1,
     &           lget
      complex*16 :: c3
      complex*16, allocatable :: a(:,:,:)
      allocate(a(2,0:maxpb,natom))
      do iatom=1,natom
        isort=is(iatom)
        i0=ip(iatom,inv_num)
        if(i0==iatom) then
          do l=0,lmpb(isort)
            a(1,l,iatom)=(0.d0,1.d0)**l
          enddo
        else if(iatom<i0) then
          do l=0,lmpb(isort)
            c3=(0.d0,1.d0)**l
            a(1,l,iatom)=c3
            a(2,l,iatom)=c3
          enddo
        else if(iatom>i0) then
          do l=0,lmpb(isort)
            c3=(0.d0,1.d0)**(l+1)
            a(1,l,iatom)=c3
            a(2,l,iatom)=c3
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
        do km1=1,n_pbmt0(isort)
          lm1=lm_pbmt(km1,isort)
          l1=lget(lm1)
          do km=1,n_pbmt0(isort)
            lm=lm_pbmt(km,isort)
            l=lget(lm)
            if(icase==1) then
              ar(km,km1,i0)=conjg(a(1,l,i0))*a(1,l1,i0)*ac(km,km1,i0)
            else if(icase==2) then
              ar(km,km1,i0)=conjg(a(1,l,i0))*ac(km,km1,i0)*a(1,l1,i0)
     &                     +conjg(a(2,l,i0))*ac(km,km1,i1)*a(2,l1,i0)
              ar(km,km1,i1)=conjg(a(1,l,i1))*ac(km,km1,i0)*a(1,l1,i1)
     &                     +conjg(a(2,l,i1))*ac(km,km1,i1)*a(2,l1,i1)
            endif
          enddo  !! over km
        enddo   !! over km1
      enddo   !! over iatom
      end
