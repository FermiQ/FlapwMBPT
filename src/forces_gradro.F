      subroutine forces_gradro(f)
	use atom_mod
	use manager_mod
      use models_mod
	use parallel_mod
	use units_mod
	use solid_mod
      use vertex_mod
      implicit none
      real*8, intent(inout) :: f(3,nsort)
      integer :: isort,m,lm,isym,mt,ispin,n,ndimv,l,irad,mm,m1,l1m1,
     &           isym1,mt1,ii,indx,icg,i,mt0
      real*8 :: dqdall,veff,const,hdif,diraddr,alf1,alf2,
     &          bet1,bet2,delt,fc
      real*8, allocatable :: work(:),work1(:),rhov(:)
      allocate(work(0:maxnrad))
      allocate(work1(0:maxnrad))
      allocate(rhov(0:maxnrad))
      const=sqrt(4.d0*pi/3.d0)
	do isort=1,nsort
        n=nrad(isort)
        hdif=1.d0/h(isort)
        ndimv=(lmpb(isort)+1)**2
        fc=4.d0*sqrt(pi)*z(isort)
        do ispin=1,nspin
          do l=0,lmpb(isort)
            do m=-l,l
              lm=l*(l+1)+m+1
              if(sym(lm,isort))then
                isym=lmsym(lm,isort)
                mt=indmt(isym,isort,ispin)
                rhov(0:n)=ro(mt:mt+n)
                if(l==0) rhov(0:n)=rhov(0:n)-ro_core(0:n,isort,ispin)
* get first-order derivative
                call dif(rhov(0),n,work1(0))
                do irad=1,n
                  diraddr=hdif/dr(irad,isort)
                  work1(irad)=work1(irad)*diraddr
     &                       -dfloat(l)*rhov(irad)/r(irad,isort)
                  work1(irad)=work1(irad)*dr(irad,isort)
     &				                   *r(irad,isort)**2
                enddo
c     *****   X - component **********************************
                mm=iabs(m)
                if(m.eq.0.or.m.eq.-1) then
                  alf1=0.d0
                else if(m.eq.1) then
                  alf1=sqrt(dfloat(l+l+1)*l*(l+1)/2.d0/(l+l-1))
                else
                  alf1=0.5d0*sqrt(dfloat(l+l+1)*(l+mm-1)*(l+mm)/(l+l-1))
                endif
                if(mm.gt.l-2) then
                  alf2=0.d0
                else if(m.eq.0) then
                  alf2=-sqrt(dfloat(l+l+1)*l*(l-1)/2.d0/(l+l-1))
                else
                 alf2=-0.5d0*sqrt(dfloat(l+l+1)*(l-mm-1)*(l-mm)/(l+l-1))
                endif
                m1=isign(mm-1,m)
                l1m1=l*(l-1)+m1+1
                if(iabs(m1).le.l-1) then
                  if(sym(l1m1,isort))then
                    isym1=lmsym(l1m1,isort)
                    mt0=indmt(isym1,isort,1)
                    mt1=indmt(isym1,isort,ispin)
                    work(0)=0.d0
                    do irad=1,n
                      veff=v_mt_h(mt0+irad)+v_mt_xc(mt1+irad)
                      if(l1m1==1) veff=veff-fc/r(irad,isort)
                      work(irad)=rhov(irad)*veff*dr(irad,isort)
     &                                          *r(irad,isort)
                    enddo
                    f(1,isort)=f(1,isort)+alf1*dqdall(h(isort),work,n)
                  endif
                endif
                m1=isign(mm+1,m)
                l1m1=l*(l-1)+m1+1
                if(iabs(m1).le.l-1) then
                  if(sym(l1m1,isort))then
                    isym1=lmsym(l1m1,isort)
                    mt0=indmt(isym1,isort,1)
                    mt1=indmt(isym1,isort,ispin)
                    work(0)=0.d0
                    do irad=1,n
                      veff=v_mt_h(mt0+irad)+v_mt_xc(mt1+irad)
                      if(l1m1==1) veff=veff-fc/r(irad,isort)
                      work(irad)=rhov(irad)*veff*dr(irad,isort)
     &                                          *r(irad,isort)
                    enddo
                    f(1,isort)=f(1,isort)+alf2*dqdall(h(isort),work,n)
                  endif
                endif
                ii=max0(lm,4)
                indx=ii*(ii-1)/2+min0(lm,4)
                do icg=indxcg(indx),indxcg(indx+1)-1
                  l1m1=jcg(icg)
                  if(l1m1.le.ndimv) then
                    if(sym(l1m1,isort)) then
                      isym1=lmsym(l1m1,isort)
                      mt0=indmt(isym1,isort,1)
                      mt1=indmt(isym1,isort,ispin)
                      work(0)=0.d0
                      do irad=1,n
                        veff=v_mt_h(mt0+irad)+v_mt_xc(mt1+irad)
                        if(l1m1==1) veff=veff-fc/r(irad,isort)
                        work(irad)=work1(irad)*veff
                      enddo
                      f(1,isort)=f(1,isort)+const
     &                          *cg(icg)*dqdall(h(isort),work,n)
                    endif
                  endif
                enddo   !!! over icg
c     *****   Y - component **********************************
                if(m.eq.0.or.m.eq.1.or.iabs(m-1).gt.l-1) then
                  bet1=0.d0
                else
                  bet1=0.5d0*sqrt(dfloat(l+l+1)*(l+m-1)*(l+m)/(l+l-1))
                  if(m.ge.0) bet1=-bet1
                endif
                if(iabs(m+1).gt.l-1) then
                  bet2=0.d0
                else if(m.eq.0) then
                  bet2=-sqrt(dfloat(l+l+1)*l*(l-1)/2.d0/(l+l-1))
                else if(m.eq.-1) then
                  bet2=sqrt(dfloat(l+l+1)*l*(l+1)/2.d0/(l+l-1))
                else
                  bet2=0.5d0*sqrt(dfloat(l+l+1)*(l-m-1)*(l-m)/(l+l-1))
                  if(m.ge.0) bet2=-bet2
                endif
                m1=-m+1
                l1m1=l*(l-1)+m1+1
                if(iabs(m1).le.l-1) then
                  if(sym(l1m1,isort))then
                    isym1=lmsym(l1m1,isort)
                    mt0=indmt(isym1,isort,1)
                    mt1=indmt(isym1,isort,ispin)
                    work(0)=0.d0
                    do irad=1,n
                      veff=v_mt_h(mt0+irad)+v_mt_xc(mt1+irad)
                      if(l1m1==1) veff=veff-fc/r(irad,isort)
                      work(irad)=rhov(irad)*veff*dr(irad,isort)
     &                                          *r(irad,isort)
                    enddo
                    f(2,isort)=f(2,isort)+bet1*dqdall(h(isort),work,n)
                  endif
                endif
                m1=-m-1
                l1m1=l*(l-1)+m1+1
                if(iabs(m1).le.l-1) then
                  if(sym(l1m1,isort))then
                    isym1=lmsym(l1m1,isort)
                    mt0=indmt(isym1,isort,1)
                    mt1=indmt(isym1,isort,ispin)
                    work(0)=0.d0
                    do irad=1,n
                      veff=v_mt_h(mt0+irad)+v_mt_xc(mt1+irad)
                      if(l1m1==1) veff=veff-fc/r(irad,isort)
                      work(irad)=rhov(irad)*veff*dr(irad,isort)
     &                                          *r(irad,isort)
                    enddo
                    f(2,isort)=f(2,isort)+bet2*dqdall(h(isort),work,n)
                  endif
                endif
                ii=max0(lm,2)
                indx=ii*(ii-1)/2+min0(lm,2)
                do icg=indxcg(indx),indxcg(indx+1)-1
                  l1m1=jcg(icg)
                  if(l1m1.le.ndimv) then
                    if(sym(l1m1,isort)) then
                      isym1=lmsym(l1m1,isort)
                      mt0=indmt(isym1,isort,1)
                      mt1=indmt(isym1,isort,ispin)
                      work(0)=0.d0
                      do irad=1,n
                        veff=v_mt_h(mt0+irad)+v_mt_xc(mt1+irad)
                        if(l1m1==1) veff=veff-fc/r(irad,isort)
                        work(irad)=work1(irad)*veff
                      enddo
                      f(2,isort)=f(2,isort)+const*cg(icg)
     &                          *dqdall(h(isort),work,n)
                    endif
                  endif
                enddo   !!! over icg
c     *****   Z - component **********************************
                if(iabs(m).gt.l-1) then
                  delt=0.d0
                else
                  delt=sqrt(dfloat(l+l+1)*(l+m)*(l-m)/(l+l-1))
                endif
                m1=m
                l1m1=l*(l-1)+m1+1
                if(iabs(m1).le.l-1) then
                  if(sym(l1m1,isort))then
                    isym1=lmsym(l1m1,isort)
                    mt0=indmt(isym1,isort,1)
                    mt1=indmt(isym1,isort,ispin)
                    work(0)=0.d0
                    do irad=1,n
                      veff=v_mt_h(mt0+irad)+v_mt_xc(mt1+irad)
                      if(l1m1==1) veff=veff-fc/r(irad,isort)
                      work(irad)=rhov(irad)*veff*dr(irad,isort)
     &                                          *r(irad,isort)
                    enddo
                    f(3,isort)=f(3,isort)+delt*dqdall(h(isort),work,n)
                  endif
                endif
                ii=max0(lm,3)
                indx=ii*(ii-1)/2+min0(lm,3)
                do icg=indxcg(indx),indxcg(indx+1)-1
                  l1m1=jcg(icg)
                  if(l1m1.le.ndimv) then
                    if(sym(l1m1,isort)) then
                      isym1=lmsym(l1m1,isort)
                      mt0=indmt(isym1,isort,1)
                      mt1=indmt(isym1,isort,ispin)
                      work(0)=0.d0
                      do irad=1,n
                        veff=v_mt_h(mt0+irad)+v_mt_xc(mt1+irad)
                        if(l1m1==1) veff=veff-fc/r(irad,isort)
                        work(irad)=work1(irad)*veff
                      enddo
                      f(3,isort)=f(3,isort)+const*cg(icg)
     &                          *dqdall(h(isort),work,n)
                    endif
                  endif
                enddo   !!! over icg
              endif  !!! over sym(lm)
            enddo   !!! over m
          enddo   !!! over l
        enddo   !!! over ispin
      enddo   !!! over isort
      deallocate(work,work1,rhov)
      if(maswrk) then
        write(iun,*)' Forces after GRADRO'
        do isort=1,nsort
          write(iun,'(i5,3(1x,f12.5))')isort,(f(i,isort),i=1,3)
        enddo
      endif
      end
