      subroutine integral_pbmt_pw_r(mt_pw,mt_pw_0,v,v0)
      use atom_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: v(3),mt_pw_0(maxlfpb,nsort),v0(3)
      real*8, intent(out) :: mt_pw(n_pbmt)
      integer :: iatom,isort,lm,l,m,ie,km,ind0,ll1,i,inv_at
      real*8 :: pi4,sqpi4,qq,pi2,phase,cc,pi4s
      real*8, allocatable :: ylm(:)
      allocate (ylm(limlpb))
      pi2=pi+pi
      pi4=4.d0*pi
      pi4s=pi4*sqrt(2.d0)
      sqpi4=sqrt(pi4)
      qq=sqrt(dot_product(v,v))
      if(qq.gt.1.d-9) then
        call sphharm(v(1)/qq,v(2)/qq,v(3)/qq,maxpb,ylm)
      else
        ylm=0.d0
        ylm(1)=1.d0/sqpi4
      endif
      do iatom=1,natom
        inv_at=ip(iatom,inv_num)
        isort=is(iatom)
        ind0=iopb(iatom)-1
        phase=pi2*dot_product(v0,tau(:,iatom))
        if(iatom==inv_at) then
          cc=pi4*cos(phase)
        else if(iatom<inv_at) then
          cc=pi4s*cos(phase)
        else
          cc=pi4s*sin(phase)
        endif
        do l=0,lmpb(isort)
          ll1=l*(l+1)+1
          do ie=1,ntle_pb(l,isort)
            i=ind_prod0(ie,l,isort)
            do m=-l,l
              lm=ll1+m
              km=ind0+ind_prod(ie,lm,isort)
              mt_pw(km)=cc*ylm(lm)*mt_pw_0(i,isort)
            enddo
          enddo  !! over ie
        enddo !! over l
      enddo   !! over iatom
      deallocate (ylm)
      end
