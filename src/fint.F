      subroutine fint(f,x,n,t,p)
c     interpolation of function f(x) in point t.
      implicit none
	integer, intent(in) :: n
      real*8, intent(in) :: f(n),x(n),t
	real*8, intent(out) :: p
	integer :: i,j
	real*8 :: s
      p=0.d0
      do i=1,n
        s=f(i)
        do j=1,i-1
          s=s*(t-x(j))/(x(i)-x(j))
        enddo
        do j=i+1,n
          s=s*(t-x(j))/(x(i)-x(j))
        enddo
        p=p+s
      enddo
      end