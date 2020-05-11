      real*8 function p_smooth_overl(ii,n,n1,l,l1,rs,st)
c --- Integral of r^ii*P(n,l)*P(n1,l1)      
      use heg_mod
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
      integer, intent(in) :: ii,n,n1,l,l1
      real*8, intent(in) :: rs,st
	integer :: nmax,i,i1,il
	real*8 :: fac2,s,ti,ti1
	real*8, allocatable :: fk(:)
	nmax=max(n,n1)
	allocate(fk(0:nmax))
	call facc(fk,nmax)
	s=0.d0
	do i1=0,n1
	  ti1=(-1)**(n1+i1)*fac2(2*l1+1)*2**i1*fk(n1)/fac2(2*(i1+l1)+1)
     &	                                /fk(i1)/fk(n1-i1)
	  do i=0,n
	    il=l+l1+2*(i+i1)
	    ti=(-1)**(n+i)*fac2(2*l+1)*2**i*fk(n)/fac2(2*(i+l)+1)/fk(i)
     &	                                         /fk(n-i)
          s=s+ti*ti1*st**(il+ii+1)/(il+ii+1.d0)/rs**il
        enddo
      enddo
	deallocate(fk)
	p_smooth_overl=s
      end