      subroutine rad_inward(e,a,b,pn,qn,l,m,r,dr,p,q,v,akap,n,nod2,q2,
     &                      key1,x_p,x_q)
	use manager_mod
	use units_mod
      IMPLICIT none
      integer, intent(in) :: l,n,key1
      real*8, intent(in) :: e,a,b,r(0:n),dr(0:n),v(0:n),akap,pn,qn
      integer, intent(inout) :: m
      real*8, intent(inout) :: x_p(0:n),x_q(0:n)
      integer, intent(out) :: nod2
      real*8, intent(out) :: p(0:n),q(0:n),q2
	integer :: i0,i,kc,ii,k,l1
      real*8 :: pp(5),qq(5),c2,r3,a3,d1,d2,d3,alf3,alf2,alf1,alf,
     &          b3,v3,r1,r2,a1,a2,b1,b2,v1,v2,vv2,xp3,xp1,xp2,
     &          ak1,am1,ak2,am2,ak3,am3,ak4,am4,xq1,xq2,xq3,
     &          adr,p_p,q_p,ppp,qqq,p_c,q_c,dqdall
      real*8, allocatable :: work(:)
	c2=clight**2
      l1=l+1
C   ASYMPTOTIC BEHAVIOR ::
      I0=n
	p(i0)=pn
	q(i0)=qn
C   4-TH ORDER RUNGE-KUTTA LOOP ::
      R3=R(I0)
      A3=DR(I0)*a
      alf3=1.d0+(e-v(i0))/c2
	if(key1.eq.1) xp3=x_p(i0)*a3
      if(irel==0) then
        B3=A3*l1/R3
        V3=A3*(v(i0)-e)
        d3=b3
      else if(irel==1) then
        B3=a3/r3
        V3=A3*(v(i0)-e+l*(l+1)/alf3/r3/r3)
        d3=b3
        a3=a3*alf3
	  if(key1.eq.1) xq3=x_q(i0)*a3/c2
      else if(irel==2) then
        B3=-a3*akap/r3
        V3=A3*(v(i0)-e)
        d3=b3
        a3=a3*alf3
	  if(key1.eq.1) xq3=x_q(i0)*a3/c2
      endif
      DO I=n-1,n-4,-1
        R1=R3
        R2=b*(EXP(a*(I+0.5D0))-1.D0) !!! R(I+1/2)
        R3=R(I)
        A1=A3
        A2=(R2+B)*a
        A3=DR(I)*a
        B1=B3
        d1=d3
        V1=V3
        alf1=alf3
        IF(I.NE.n-1) Vv2=(9.D0*(v(I)+v(I+1))-v(I-1)-v(I+2))*6.25D-2
        IF(I.EQ.n-1) vV2=(v(I)+v(I+1))/2
	  if(key1.eq.1) then
	    xp1=xp3
          IF(I.NE.n-1) xp2=A2*
     &		(9.D0*(x_p(I)+x_p(I+1))-x_p(I-1)-x_p(I+2))*6.25D-2
          IF(I.EQ.n-1) xp2=A2*(x_p(I)+x_p(I+1))*0.5D0
	    xp3=x_p(i)*a3
          if(irel/=0) then
	      xq1=xq3
            IF(I.NE.n-1) xq2=a2/c2*
     &		(9.D0*(x_q(I)+x_q(I+1))-x_q(I-1)-x_q(I+2))*6.25D-2
            IF(I.EQ.n-1) xq2=a2/c2*(x_q(I)+x_q(I+1))*0.5D0
	      xq3=x_q(i)*a3/c2
          endif
	  endif
        if(irel==0) then
          B2=A2/R2*L1
          B3=A3/R3*L1
          d2=a2*l1/r2
          d3=a3*l1/r3
          V2=A2*(vv2-e)
          V3=A3*(v(I)-e)
        else if(irel==1) then
          alf2=1.d0+(e-vv2)/c2
          alf3=1.d0+(e-v(i))/c2
          B2=a2/r2
          B3=a3/r3
          d2=a2/r2
          d3=a3/r3
          V2=A2*(vv2-e+l*(l+1)/alf2/r2/r2)
          V3=A3*(v(I)-e+l*(l+1)/alf3/r3/r3)
          a2=a2*alf2
          a3=a3*alf3
        else if(irel==2) then
          alf2=1.d0+(e-vv2)/c2
          alf3=1.d0+(e-v(i))/c2
          B2=-a2*akap/r2
          B3=-a3*akap/r3
          d2=b2
          d3=b3
          V2=A2*(vv2-e)
          V3=A3*(v(I)-e)
          a2=a2*alf2
          a3=a3*alf3
        endif
C   **************
        aK1= A1*Q(I+1)+B1*P(I+1)
        aM1=V1*P(I+1)-d1*Q(I+1)
	  if(key1.eq.1) am1=am1+xp1
	  if(key1==1.and.irel/=0) ak1=ak1-xq1
C   **************
        aK2= A2*(Q(I+1)-aM1*0.5D0)+B2*(P(I+1)-aK1*0.5D0)
        aM2=V2*(P(I+1)-aK1*0.5D0)-d2*(Q(I+1)-aM1*0.5D0)
	  if(key1.eq.1) am2=am2+xp2
	  if(key1==1.and.irel/=0) ak2=ak2-xq2
C   **************
        aK3= A2*(Q(I+1)-aM2*0.5D0)+B2*(P(I+1)-aK2*0.5D0)
        aM3= V2*(P(I+1)-aK2*0.5D0)-d2*(Q(I+1)-aM2*0.5D0)
	  if(key1.eq.1) am3=am3+xp2
	  if(key1==1.and.irel/=0) ak3=ak3-xq2
C   **************
        aK4= A3*(Q(I+1)-aM3)+B3*(P(I+1)-aK3)
        aM4= V3*(P(I+1)-aK3)-d3*(Q(I+1)-aM3)
	  if(key1.eq.1) am4=am4+xp3
	  if(key1==1.and.irel/=0) ak4=ak4-xq3
C   **************
        P(I)=P(I+1)-(aK1+aK4+2.D0*(aK2+aK3))/6.D0
        Q(I)=Q(I+1)-(aM1+aM4+2.D0*(aM2+aM3))/6.D0
      enddo


	do i=1,5
	  ii=n-i+1
	  adr=a*dr(ii)
        if(irel==0) then
	    pp(i)=l1*p(ii)/r(ii)+q(ii)
	    qq(i)=(v(ii)-e)*p(ii)-l1*q(ii)/r(ii)
        else if(irel==1) then
          alf=1.d0+(e-v(ii))/c2
	    pp(i)=alf*q(ii)+p(ii)/r(ii)
	    qq(i)=(v(ii)-e+l*(l+1)/alf/r(ii)**2)*p(ii)-q(ii)/r(ii)
	    if(key1.eq.1) pp(i)=pp(i)-x_q(ii)/c2
        else if(irel==2) then
          alf=1.d0+(e-v(ii))/c2
	    pp(i)=alf*q(ii)-akap/r(ii)*p(ii)
	    qq(i)=(v(ii)-e)*p(ii)+akap*q(ii)/r(ii)
	    if(key1.eq.1) pp(i)=pp(i)-x_q(ii)/c2
        endif
	  pp(i)=adr*pp(i)
	  if(key1.eq.1) qq(i)=qq(i)+x_p(ii)
	  qq(i)=adr*qq(i)
	enddo

	do i=n-4,m+1,-1
	  p_p=p(i)-(1901.d0*pp(5)-2774.d0*pp(4)+2616.d0*pp(3)
     &	       -1274.d0*pp(2)+251.d0*pp(1))/720.d0
	  q_p=q(i)-(1901.d0*qq(5)-2774.d0*qq(4)+2616.d0*qq(3)
     &	       -1274.d0*qq(2)+251.d0*qq(1))/720.d0
	  adr=a*dr(i-1)
        if(irel==0) then
	    ppp=l1*p_p/r(i-1)+q_p
	    qqq=(v(i-1)-e)*p_p-l1*q_p/r(i-1)
        else if(irel==1) then
          alf=1.d0+(e-v(i-1))/c2
	    ppp=alf*q_p+p_p/r(i-1)
	    qqq=(v(i-1)-e+l*(l+1)/alf/r(i-1)**2)*p_p-q_p/r(i-1)
	    if(key1.eq.1) ppp=ppp-x_q(i-1)/c2
        else if(irel==2) then
          alf=1.d0+(e-v(i-1))/c2
	    ppp=alf*q_p-akap/r(i-1)*p_p
	    qqq=(v(i-1)-e)*p_p+akap*q_p/r(i-1)
	    if(key1.eq.1) ppp=ppp-x_q(i-1)/c2
        endif
	  ppp=adr*ppp
	  if(key1.eq.1) qqq=qqq+x_p(i-1)
	  qqq=adr*qqq
	  p_c=p(i)-(251.d0*ppp+646.d0*pp(5)-264.d0*pp(4)
     &	        +106.d0*pp(3)-19.d0*pp(2))/720.d0
	  p(i-1)=(475.d0*p_c+27.d0*p_p)/502.d0
	  q_c=q(i)-(251.d0*qqq+646.d0*qq(5)-264.d0*qq(4)
     &	        +106.d0*qq(3)-19.d0*qq(2))/720.d0
	  q(i-1)=(475.d0*q_c+27.d0*q_p)/502.d0
	  do k=1,4
	    pp(k)=pp(k+1)
	    qq(k)=qq(k+1)
	  enddo
	  adr=a*dr(i-1)
        if(irel==0) then
	    pp(5)=l1*p(i-1)/r(i-1)+q(i-1)
	    qq(5)=(v(i-1)-e)*p(i-1)-l1*q(i-1)/r(i-1)
        else if(irel==1) then
          alf=1.d0+(e-v(i-1))/c2
	    pp(5)=alf*q(i-1)+p(i-1)/r(i-1)
	    qq(5)=(v(i-1)-e+l*(l+1)/alf/r(i-1)**2)*p(i-1)-q(i-1)/r(i-1)
	    if(key1.eq.1) pp(5)=pp(5)-x_q(i-1)/c2
        else if(irel==2) then
          alf=1.d0+(e-v(i-1))/c2
	    pp(5)=alf*q(i-1)-akap/r(i-1)*p(i-1)
	    qq(5)=(v(i-1)-e)*p(i-1)+akap*q(i-1)/r(i-1)
	    if(key1.eq.1) pp(5)=pp(5)-x_q(i-1)/c2
        endif
	  pp(5)=adr*pp(5)
	  if(key1.eq.1) qq(5)=qq(5)+x_p(i-1)
	  qq(5)=adr*qq(5)
	  if(mod(i-1,2).eq.0) then
		  if(p(i-1)*pp(5).gt.0.d0) then
	      kc=i-1
            goto 2
	    endif
	  endif
	enddo
	kc=m
2	continue
	m=kc


	nod2=0
	do i=m,n-5
	  if(p(i)*p(i+1).lt.0.d0) nod2=nod2+1
	enddo
      allocate(work(m:n))
	do i=m,n
	  work(i)=p(i)**2*dr(i)
      enddo
	if(irel.ge.1) then
	  do i=m,n
	    work(i)=work(i)+q(i)**2*dr(i)/c2
        enddo
	endif
      q2=dqdall(a,work(m),n-m)
      deallocate(work)
      END