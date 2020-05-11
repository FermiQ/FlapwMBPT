      subroutine mathpaw(nsitmax,radius,rofi,roatom)
***************************************************************
*  Perform a summation of the pseudo charge tails over lattice*
***************************************************************
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: nsitmax
      real*8, intent(in) :: roatom(0:nrmax,nsort),rofi(0:nrmax,nsort),
     & 	                    radius
      integer :: isort,nvcmax,ierr,mt,irad,i,i0,isite,jatom,jsort,irad1,
     &	         irad2,num,ilmax,ii,l,m,lm,isym,lp,isymp,ind_site
      real*8 :: sqpi4,rovtt,rmd,rpd,hout1,xpnt
      integer, allocatable :: iref(:,:),nsite(:),nd(:),n0(:)
      real*8,allocatable :: work(:,:),answ(:),pol(:),rnew(:),delta(:,:),
     &	                    rovt1(:),wd(:,:,:),rho_sort(:)
      allocate(work(0:maxnrad,maxnsym),answ(maxnsym),pol(0:maxpb),
     &	       rnew(0:maxnrad),delta(nsitmax,nsort),rovt1(0:maxnrad),
     &         wd(nsitmax,maxnsym,nsort))
      allocate(iref(nsitmax,nsort),nsite(nsort))
      nvcmax=30000
      sqpi4=2.d0*sqrt(pi)
      call setsites(delta,iref,nsite,nvcmax,wd,nsitmax,ierr,radius)
      if(ierr.gt.0) then
        if(maswrk) write(iun,*)' Exit from Setsites in Mathpaw : Err'
        call ending
      endif
c *** calculate tail contribution to the pseudo density ***
      allocate(nd(nproc))
      allocate(n0(nproc))
      do isort=1,nsort
        mt=indmt(1,isort,1)
        call dcopy(nrad(isort)+1,roatom(0,isort),1,ro(mt),1)
        do irad=1,nrad(isort)
          do i=nr(isort),0,-1
            if(rofi(i,isort).lt.r(irad,isort)) then
              i0=max(0,i-1)
              goto 50
            endif
          enddo
50        call fint(roatom(i0,isort),rofi(i0,isort),4,r(irad,isort),
     &              ro(mt+irad))
        enddo
        call size_shift_par(nsite(isort),nproc,nd,n0)
        allocate(rho_sort(0:nrad(isort)))
        rho_sort=0.d0
        do ind_site=1,nd(me+1)
          isite=n0(me+1)+ind_site
          jatom=iref(isite,isort)
          jsort=is(jatom)
          if(delta(isite,isort).eq.0.d0) goto 80
          do i=1,nr(jsort)
            if(rofi(i,jsort).gt.delta(isite,isort)) then
              i0=i-2
              goto 5
            endif
          enddo
5         call fint(roatom(i0,jsort),rofi(i0,jsort),4,
     &              delta(isite,isort),rovtt)
          mt=indmt(1,isort,1)
          rho_sort(0)=rho_sort(0)+sqpi4*rovtt*wd(isite,1,isort)
          do irad=1,nrad(isort)
c  make a grid for taking the integrals
            rmd=delta(isite,isort)-r(irad,isort)
            rpd=delta(isite,isort)+r(irad,isort)
            do i=1,nr(jsort)
              if(rofi(i,jsort).gt.rmd) then
                irad1=i-1
                goto 6
              endif
            enddo
6           do i=1,nr(jsort)
              if(rofi(i,jsort).gt.rpd) then
                irad2=i
                goto 7
              endif
            enddo
7           num=irad2-irad1+60
c            ilmax=lmpb(isort)
	        ilmax=0
            hout1=(rpd-rmd)/num
c  make the integrand
            do i=0,num
              rnew(i)=rmd+i*hout1
              do ii=1,nr(jsort)
                if(rofi(ii,jsort).gt.rnew(i)) then
                  i0=ii-2
                  goto 8
                endif
              enddo
8             call fint(roatom(i0,jsort),rofi(i0,jsort),4,rnew(i),
     &                  rovt1(i))
              xpnt=0.5d0*(r(irad,isort)**2+delta(isite,isort)**2-
     &             rnew(i)**2)/(r(irad,isort)*delta(isite,isort))
              call legpol(pol,xpnt,ilmax)
              do l=0,ilmax
                do m=-l,l
                  lm=l*(l+1)+m+1
                  if(sym(lm,isort)) then
                    isym=lmsym(lm,isort)
                    work(i,isym)=rnew(i)*pol(l)
                  endif
                enddo
              enddo
            enddo     !!! over i
            lp=-1
            do l=0,ilmax
              do m=-l,l
                lm=l*(l+1)+m+1
                if(sym(lm,isort)) then
                  isym=lmsym(lm,isort)
                  if(lp.ne.l) then
                    do i=0,num
                      work(i,isym)=work(i,isym)*rovt1(i)
                    enddo
                    answ(isym)=0.d0
                    do ii=0,num-1
                     answ(isym)=answ(isym)+work(ii,isym)+work(ii+1,isym)
                    enddo
                    answ(isym)=answ(isym)*hout1/2.d0
                    lp=l
                    isymp=isym
                  else
                    answ(isym)=answ(isymp)
                  endif
                endif
              enddo  !!! over m
            enddo     !!! over l
            do l=0,ilmax
              do m=-l,l
                lm=l*(l+1)+m+1
                if(sym(lm,isort)) then
                  isym=lmsym(lm,isort)
                  rho_sort(irad)=rho_sort(irad)+2.d0*pi/r(irad,isort)
     &				/sqpi4/delta(isite,isort)*answ(isym)*
     &                     wd(isite,isym,isort)
                endif
              enddo
            enddo
          enddo   !!! over irad
80        continue
        enddo   !!! over isite
        if(nproc/=1) then
#ifdef MPI
          call DGOP(rho_sort,nrad(isort)+1,'  +',MPI_COMM_WORLD)
#endif
        endif
        mt=indmt(1,isort,1)
        do irad=0,nrad(isort)
          ro(mt+irad)=ro(mt+irad)+rho_sort(irad)
        enddo
        deallocate(rho_sort)
      enddo   !!! over isort
      deallocate(work,answ,pol,rnew,delta,rovt1,wd,iref,nsite,nd,n0)
      end
