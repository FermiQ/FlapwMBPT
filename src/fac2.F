      real*8 function fac2(l)
	integer :: l,i,ll
      if(l.eq.-1) then
        fac2=1
        return
      endif
      if(l.eq.1) then
        fac2=1
        return
      endif
      fac2=1
      ll=(l-1)/2
      do i=1,ll
        fac2=fac2*(2*i+1)
	enddo
      end
      
      subroutine facc(f,n)
      implicit none
	integer, intent(in) :: n
      real*8, intent(out) :: f(0:n)
	integer :: i
      f(0)=1.d0
      do i=1,n
        f(i)=f(i-1)*i
      enddo
      end
      
      real*8 function fact12(n,m)   !!!   n!/m!
c***************************************************
      implicit none
	integer, intent(in) :: n,m
	integer :: i
      if(n.lt.m)then
        fact12=1.d0/dble(n+1)
        do i=n+2,m
          fact12=fact12/dble(i)
        enddo
      elseif(n.eq.m)then
        fact12=1.d0
      elseif(n.gt.m)then
        fact12=dble(m+1)
        do i=m+2,n
          fact12=fact12*dble(i)
        enddo
      endif
      end
      
      real*8 function gamma_f(x)
c***** Gamma function of integer or half-integer ***************
      use units_mod
      implicit none
	real*8, intent(in) :: x
	integer :: i,n,j
      real*8 :: f
      i=idnint(x+x)
      if(mod(i,2)==0) then
        n=i/2
        f=1.d0
        do j=2,n-1
          f=f*j
        enddo
      else
        n=i
        f=1.d0
        do j=2,n-1
          f=f*j
        enddo
        do j=2,(n-1)/2
          f=f/j
        enddo
        f=f*sqrt(pi)/2.d0**(n-1)
      endif
      gamma_f=f
      end
      
      real*8 function hyper_1f2(a1,b1,b2,z)
c***** Gamma function of integer or half-integer ***************
      use units_mod
      implicit none
	real*8, intent(in) :: a1,b1,b2,z
	integer :: k
      real*8 :: fk,s
      fk=a1/b1/b2
      s=fk
      k=0
      do while(abs(fk)>1.d-30)
        k=k+1
        fk=fk*(a1+k)*z/(b1+k)/(b2+k)
        s=s+fk
      enddo
      hyper_1f2=s
      end
      
      real*8 function hyper_1f2_reg(a1,b1,b2,z)
c***** Gamma function of integer or half-integer ***************
      use units_mod
      implicit none
	real*8, intent(in) :: a1,b1,b2,z
      real*8 :: gamma_f,hyper_1f2
      hyper_1f2_reg=hyper_1f2(a1,b1,b2,z)/gamma_f(b1)/gamma_f(b2)
      end