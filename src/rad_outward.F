       subroutine rad_outward(e,a,l,ii,kc,r,dr,p,q,v,akap,n,nod1,q1,ai,
     &                        bi,z,key,p0,q0,p2,key1,x_p,x_q)
	use manager_mod
	use units_mod
      IMPLICIT none
      integer, intent(in) :: l,n,ii,key,key1
      real*8, intent(in) :: e,a,r(0:n),dr(0:n),v(0:n),akap,ai(0:3),
     &                      bi(0:3),z,p0(0:n),q0(0:n),p2(0:n)
      integer, intent(inout) :: kc
      real*8, intent(inout) :: x_p(0:n),x_q(0:n)
      integer, intent(out) :: nod1
      real*8, intent(out) :: p(0:n),q(0:n),q1
	integer :: i,k,l1
      real*8 :: pp(5),qq(5),c2,alf,ppp,adr,p_p,q_p,p_c,q_c,dqdall,qqq
      real*8, allocatable :: work(:)
      l1=l+1
	c2=clight**2
C   ASSYMPTOTIC R ----> 0   BEHAVIOR ::
      call pq_r0_asym(5,l,ii,r,ai,bi,p,q,z)
c -----------------------------------------------------------------      
      do i=1,5
	  adr=a*dr(i)
        if(irel==0) then
          pp(i)=l1*p(i)/r(i)+q(i)
	    qq(i)=(v(i)-e)*p(i)-l1*q(i)/r(i)
          if(key==1) qq(i)=qq(i)-p0(i)
          if(key==2) qq(i)=qq(i)-2*p0(i)
        else if(irel==1) then
          alf=1.d0+(e-v(i))/c2
	    pp(i)=alf*q(i)+p(i)/r(i)
	    if(key1.eq.1) pp(i)=pp(i)-x_q(i)/c2
	    qq(i)=(v(i)-e+l*(l+1)/alf/r(i)**2)*p(i)-q(i)/r(i)
          if(key==1) then
            pp(i)=pp(i)+q0(i)/c2
            qq(i)=qq(i)-(1.d0+l*(l+1)/alf**2/r(i)**2/c2)*p0(i)
          endif
          if(key==2) then
            pp(i)=pp(i)+2*q0(i)/c2
            qq(i)=qq(i)-2*(1.d0+l*(l+1)/alf**2/r(i)**2/c2)*p0(i)
     &           +2*l*(l+1)/alf**3/r(i)**2/c2/c2*p2(i)
          endif
        else if(irel==2) then
          alf=1.d0+(e-v(i))/c2
	    pp(i)=alf*q(i)-akap/r(i)*p(i)
	    if(key1.eq.1) pp(i)=pp(i)-x_q(i)/c2
	    qq(i)=(v(i)-e)*p(i)+akap*q(i)/r(i)
          if(key==1) then
            pp(i)=pp(i)+q0(i)/c2
            qq(i)=qq(i)-p0(i)
          endif
          if(key==2) then
            pp(i)=pp(i)+2*q0(i)/c2
            qq(i)=qq(i)-2*p0(i)
          endif
        endif
	  pp(i)=adr*pp(i)
	  if(key1.eq.1) qq(i)=qq(i)+x_p(i)
	  qq(i)=adr*qq(i)
	enddo

	do i=5,kc-1
	  p_p=p(i)+(1901.d0*pp(5)-2774.d0*pp(4)+2616.d0*pp(3)
     &	       -1274.d0*pp(2)+251.d0*pp(1))/720.d0
	  q_p=q(i)+(1901.d0*qq(5)-2774.d0*qq(4)+2616.d0*qq(3)
     &	       -1274.d0*qq(2)+251.d0*qq(1))/720.d0
	  adr=a*dr(i+1)
        if(irel==0) then
	    ppp=l1*p_p/r(i+1)+q_p
	    qqq=(v(i+1)-e)*p_p-l1*q_p/r(i+1)
          if(key==1) qqq=qqq-p0(i+1)
          if(key==2) qqq=qqq-2*p0(i+1)
        else if(irel==1) then
          alf=1.d0+(e-v(i+1))/c2
	    ppp=alf*q_p+p_p/r(i+1)
	    if(key1.eq.1) ppp=ppp-x_q(i+1)/c2
	    qqq=(v(i+1)-e+l*(l+1)/alf/r(i+1)**2)*p_p-q_p/r(i+1)
          if(key==1) then
            ppp=ppp+q0(i+1)/c2
            qqq=qqq-(1.d0+l*(l+1)/alf**2/r(i+1)**2/c2)*p0(i+1)
          endif
          if(key==2) then
            ppp=ppp+2*q0(i+1)/c2
            qqq=qqq-2*(1.d0+l*(l+1)/alf**2/r(i+1)**2/c2)*p0(i+1)
     &             +2*l*(l+1)/alf**3/r(i+1)**2/c2/c2*p2(i+1)
          endif
        else if(irel==2) then
          alf=1.d0+(e-v(i+1))/c2
	    ppp=alf*q_p-akap/r(i+1)*p_p
	    if(key1.eq.1) ppp=ppp-x_q(i+1)/c2
	    qqq=(v(i+1)-e)*p_p+akap*q_p/r(i+1)
          if(key==1) then
            ppp=ppp+q0(i+1)/c2
            qqq=qqq-p0(i+1)
          endif
          if(key==2) then
            ppp=ppp+2*q0(i+1)/c2
            qqq=qqq-2*p0(i+1)
          endif
        endif
	  ppp=adr*ppp
	  if(key1.eq.1) qqq=qqq+x_p(i+1)
	  qqq=adr*qqq
	  p_c=p(i)+(251.d0*ppp+646.d0*pp(5)-264.d0*pp(4)
     &	       +106.d0*pp(3)-19.d0*pp(2))/720.d0
	  p(i+1)=(475.d0*p_c+27.d0*p_p)/502.d0
	  q_c=q(i)+(251.d0*qqq+646.d0*qq(5)-264.d0*qq(4)
     &	       +106.d0*qq(3)-19.d0*qq(2))/720.d0
	  q(i+1)=(475.d0*q_c+27.d0*q_p)/502.d0
	  do k=1,4
	    pp(k)=pp(k+1)
	    qq(k)=qq(k+1)
	  enddo
	  adr=a*dr(i+1)
        if(irel==0) then
          pp(5)=l1*p(i+1)/r(i+1)+q(i+1)
	    qq(5)=(v(i+1)-e)*p(i+1)-l1*q(i+1)/r(i+1)
          if(key==1) qq(5)=qq(5)-p0(i+1)
          if(key==2) qq(5)=qq(5)-2*p0(i+1)
        else if(irel==1) then
          alf=1.d0+(e-v(i+1))/c2
	    pp(5)=alf*q(i+1)+p(i+1)/r(i+1)
	    if(key1.eq.1) pp(5)=pp(5)-x_q(i+1)/c2
	    qq(5)=(v(i+1)-e+l*(l+1)/alf/r(i+1)**2)*p(i+1)-q(i+1)/r(i+1)
          if(key==1) then
            pp(5)=pp(5)+q0(i+1)/c2
            qq(5)=qq(5)-(1.d0+l*(l+1)/alf**2/r(i+1)**2/c2)*p0(i+1)
          endif
          if(key==2) then
            pp(5)=pp(5)+2*q0(i+1)/c2
            qq(5)=qq(5)-2*(1.d0+l*(l+1)/alf**2/r(i+1)**2/c2)*p0(i+1)
     &                 +2*l*(l+1)/alf**3/r(i+1)**2/c2/c2*p2(i+1)
          endif
        else if(irel==2) then
          alf=1.d0+(e-v(i+1))/c2
	    pp(5)=alf*q(i+1)-akap/r(i+1)*p(i+1)
	    if(key1.eq.1) pp(5)=pp(5)-x_q(i+1)/c2
	    qq(5)=(v(i+1)-e)*p(i+1)+akap*q(i+1)/r(i+1)
          if(key==1) then
            pp(5)=pp(5)+q0(i+1)/c2
            qq(5)=qq(5)-p0(i+1)
          endif
          if(key==2) then
            pp(5)=pp(5)+2*q0(i+1)/c2
            qq(5)=qq(5)-2*p0(i+1)
          endif
        endif
	  pp(5)=adr*pp(5)
	  if(key1.eq.1) qq(5)=qq(5)+x_p(i+1)
	  qq(5)=adr*qq(5)
	enddo

	nod1=0
	do i=0,kc-1
	  if(p(i)*p(i+1).lt.0.d0) nod1=nod1+1
	enddo
      allocate(work(0:kc))
	do i=0,kc
	  work(i)=p(i)**2*dr(i)
      enddo
	if(irel.ge.1) then
	  do i=0,kc
	    work(i)=work(i)+q(i)**2*dr(i)/c2
        enddo
	endif
      q1=dqdall(a,work(0),kc)
      deallocate(work)
      END