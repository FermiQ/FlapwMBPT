      subroutine dp_dq(e,z,l,r,p,q,v,akap,p0,q0,key2,n,p1,dp,dq)
c        Radial derivatives
c	  key2 = 1  For energy derivative
c	  key2 = 2  For second energy derivative
	  use manager_mod
	  use units_mod
      IMPLICIT none
	  integer, intent(in) :: n,key2,l
      real*8, intent(in) :: e,z,r(0:n),p(0:n),q(0:n),v(0:n),akap,
     &                      p0(0:n),q0(0:n),p1(0:n)
      real*8, intent(out) :: dp(0:n),dq(0:n)
	  integer :: l1,i
      real*8 :: c2,z2,alf,ev,fll,r2
      L1=L+1
	  c2=clight**2
	  z2=z+z
      dp(0)=0.d0
      dq(0)=0.d0
	  do i=1,n
        ev=e+z2/r(i)-v(i)
	    r2=r(i)**2
	    if(irel.eq.0) then
		  dp(i)=l1*p(i)/r(i)+q(i)
	      dq(i)=-ev*p(i)-l1*q(i)/r(i)
		  if(key2.ge.1) dq(i)=dq(i)-key2*p0(i)
	    else if(irel.eq.1) then
	      alf=1.d0+ev/c2
          fll=akap/r(i)**2/alf
		  dp(i)=p(i)/r(i)+alf*q(i)
	      dq(i)=(fll-ev)*p(i)-q(i)/r(i)
		  if(key2.ge.1) then
	        dp(i)=dp(i)+key2*q0(i)/c2
		    dq(i)=dq(i)-key2*(1.d0+fll/alf/c2)*p0(i)
		    if(key2.eq.2) dq(i)=dq(i)+2.d0*fll/alf**2/c2**2*p1(i)
          endif
	    else if(irel.eq.2) then
	      alf=1.d0+ev/c2
		  dp(i)=-akap*p(i)/r(i)+alf*q(i)
	      dq(i)=-ev*p(i)+akap*q(i)/r(i)
		  if(key2.ge.1) then
	        dp(i)=dp(i)+key2*q0(i)/c2
		    dq(i)=dq(i)-key2*p0(i)
          endif
	    endif
	  enddo
      END
