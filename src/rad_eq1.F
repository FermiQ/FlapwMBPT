      subroutine rad_eq1(e,tol,l,i,ptn,v0,fifi,pnew,qnew,work,isort,
     &                   nrad0,nod1,phi,dphi,dny1)
c	  The outward solution for given T=1/2-1/pi*arctan(D)+N
c	  Input:
c		ptn = potential function (see above)
c	  Output:
c		pnew = P_new(r)
c		qnew = Q_new(r)
	  use atom_mod
	  use manager_mod
	  use units_mod
      implicit none
      integer, intent(in) :: l,i,nrad0,isort
      real*8, intent(in) :: tol,ptn,v0(0:nrad0)
      real*8, intent(inout) :: pnew(0:nrad0),qnew(0:nrad0),
     &                         work(0:nrad0)
      integer, intent(out) :: nod1
      real*8, intent(out) :: e,fifi,dny1,phi,dphi
	  integer :: icase,it,nod,nn
      real*8 :: hh,dff,tny,tny0,e0,e1,f0,f1,f00
      real*8, allocatable :: p1(:),q1(:)
	  allocate(p1(0:nrad0),q1(0:nrad0))
	  nn=int(ptn)
	  nod=nn-l-1
	  tny0=ptn-dfloat(nn)+dfloat(nod)
	  icase=1
	  hh=0.5d0
	  it=0
	  dff=1.d0
	  do while(it<=200.and.dff>tol)
        call rad_eq0(e,z(isort),l,i,nod1,v0,fifi,r(0,isort),nrad(isort),
     &	             pnew,qnew,phi,dphi,dny1,dr(0,isort),0,h(isort),
     &               p1,q1,p1,q1,work)
c
c       the beginning of the iteration loop.
c
	  if(icase==1) then
	    e0=e
		tny=0.5d0+dfloat(nod1)-atan(dny1)/pi
	    f0=tny-tny0
	    icase=2
	    if(abs(f0)<tol) then
	      goto 5
	    else if(f0<0.d0) then
	      hh=abs(hh)
	    else
	      hh=-abs(hh)
	    endif
	    e=e+hh
	  else if(icase==2) then
		tny=0.5d0+dfloat(nod1)-atan(dny1)/pi
	    f00=tny-tny0
	    if(f0*f00<0.d0) then
	      icase=3
	      f1=f00
	      e1=e
	      e=0.5d0*(e0+e1)
	    else
	      f0=f00
	      e0=e
	      e=e+hh
	    endif
	  else if(icase==3) then
		tny=0.5d0+dfloat(nod1)-atan(dny1)/pi
	    f00=tny-tny0
	    if(f0*f00<0.d0) then
	      f1=f00
	      e1=e
	    else
	      f0=f00
	      e0=e
	    endif
	    e=0.5d0*(e0+e1)
	    dff=abs(e1-e0)
	  endif
	  it=it+1
	  enddo
5	  deallocate(p1,q1)
      end
      
      
      subroutine rad_eq4(e,tol,l,i,nodes,v0,fifi,pnew,qnew,
     &                   work,isort,nrad0,node,phi,dphi,dny1)
c	  The outward solution for the bound solution P(S)=0
c	  Input:
c	  Output:
c		pnew = P_new(r)
c		qnew = Q_new(r)
	  use atom_mod
	  use manager_mod
	  use units_mod
      implicit none
      integer, intent(in) :: l,i,nrad0,isort,nodes
      real*8, intent(in) :: tol,v0(0:nrad0)
      real*8, intent(inout) :: pnew(0:nrad0),qnew(0:nrad0),work(0:nrad0)
      integer, intent(out) :: node
      real*8, intent(out) :: e,fifi,phi,dphi,dny1
	  integer :: icase,it
      real*8 :: hh,dff,e0,e1,f0,f1,f
      real*8, allocatable :: p1(:),q1(:)
	  allocate(p1(0:nrad0),q1(0:nrad0))
	  icase=1
	  hh=0.05d0
	  it=0
	  dff=1.d0
	  e0=-1.d15
	  e1=1.d15
	  do while(it<=100.and.dff>tol)
        call rad_eq0(e,z(isort),l,i,node,v0,fifi,r(0,isort),nrad(isort),
     &	             pnew,qnew,phi,dphi,dny1,dr(0,isort),0,h(isort),
     &               p1,q1,p1,q1,work)
	  f=pnew(nrad0)
	  if(icase==1) then
	    if(abs(f)<tol) then
	      goto 5
	    else if(node==nodes) then
	      if(e>e0) then
	        e0=e
	        f0=f
	      endif
	      e=e+hh
	    else if(node<nodes) then
	      e=e+hh
	    else if(node>nodes) then
	      if(e<e1) then
	        e1=e
	        f1=f
	      endif
	      e=e-hh
	    endif
	    if(e0>-1.d14.and.e1<1.d14) then
	      icase=2
	      e=0.5d0*(e0+e1)
	    endif
	  else if(icase==2) then
	    if(abs(f)<tol) then
	      goto 5
	    else if(node==nodes.and.f*f0>=0.d0) then
	        e0=e
	    else if(node>nodes.or.f*f1<0.d0) then
	        e1=e
	    endif
	    e=0.5d0*(e0+e1)
	  endif
	  it=it+1
	  enddo
5     deallocate(p1,q1)
      end
