      subroutine geom_factor_finder(n,a,q)
c     a = Full range / first step
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: a
      real*8, intent(out) :: q
      integer :: m
	  real*8 :: q1,q2,q3,f1,f2,f3,f_rad_fun
      q1=1.d-8
      q2=3.d0
      f1=f_rad_fun(q1,a,n)
      f2=f_rad_fun(q2,a,n)
      m=0
1     q3=(q1+q2)/2
      f3=f_rad_fun(q3,a,n)
      if(f3<0.d0) then
        q1=q3
        f1=f3
      else
        q2=q3
        f2=f3
      endif
      m=m+1
      if(f2-f1>1.d-12.and.m<60) goto 1
      q=1.d0+(q1+q2)/2
      if(maswrk.and.m==60) write(iun,*)
     &        'GEOM_FACTOR_FINDER: f2-f1 is ',f2-f1
      end



      subroutine geom_factor_finder_rad(n,a,q)
c     a = Full range / first step
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: a
      real*8, intent(out) :: q
      integer :: m
	  real*8 :: q1,q2,q3,f1,f2,f3,f_rad_fun
      q1=1.d-8
      q2=0.1d0
      f1=f_rad_fun(q1,a,n)
      f2=f_rad_fun(q2,a,n)
      m=0
1     q3=(q1+q2)/2
      f3=f_rad_fun(q3,a,n)
      if(f3<0.d0) then
        q1=q3
        f1=f3
      else
        q2=q3
        f2=f3
      endif
      m=m+1
      if(f2-f1>1.d-11.and.m<60) goto 1
      q=1.d0+(q1+q2)/2
      if(maswrk.and.m==60) write(iun,*)
     &        'GEOM_FACTOR_FINDER: f2-f1 is ',f2-f1
      end

      real*8 function f_rad_fun(x,a,n)
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: x,a
      integer :: k
      real*8 :: f,f0
      f0=-a/n
      f=1.d0
      do k=1,10000
        f0=f0+f
        f=f*x*(n-k)/dfloat(k+1)
        if(f<1.d-30) exit
      enddo
      f_rad_fun=f0
      end
