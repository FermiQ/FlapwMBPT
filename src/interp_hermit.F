      subroutine interp_hermit(x0,x1,y0,y1,dy0,dy1,x,y)
c     interpolation of function f(x) in point t.
      implicit none
      real*8, intent(in) :: x0,x1,y0,y1,dy0,dy1,x
	real*8, intent(out) :: y
	real*8 :: h,t,t2,t3
      h=x1-x0
      t=(x-x0)/h
      t2=t*t
      t3=t2*t
      y=y0*(2*t3-3*t2+1)+y1*(-2*t3+3*t2)+h*dy0*(t3-2*t2+t)+h*dy1*(t3-t2)
      end