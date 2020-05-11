      subroutine spline_inhmg_c(mesh,y,a,b,c,d,n,key,dy_0,dy_n)
c	key = 0 - Free boundary conditions (y''(0)=y''(n)=0)
c	key = 1 - y'(0) and y'(n) are given
c	key = 2 - y''(0)=0 and y'(n) is given
c	key = 3 - y'(0) is given and y''(n)=0
      implicit none
	integer, intent(in) :: n,key
	real*8, intent(in) :: mesh(0:n)
	complex*16, intent(in) :: y(0:n),dy_0,dy_n
	complex*16, intent(out) :: a(n),b(n),c(n),d(n)
      real*8 :: d0,dn
	real*8, allocatable :: f(:),aa(:,:)
	allocate(f(0:n))
      allocate(aa(n,4))
c --------- Real part ---------------------------------------
      f=real(y)
      d0=real(dy_0)
      dn=real(dy_n)
      call spline_inhmg(mesh,f,aa(1,1),aa(1,2),aa(1,3),aa(1,4),n,key,d0,
     &                  dn)
      a=aa(:,1)
      b=aa(:,2)
      c=aa(:,3)
      d=aa(:,4)
c --------- Imaginary part ---------------------------------------
      f=imag(y)
      d0=imag(dy_0)
      dn=imag(dy_n)
      call spline_inhmg(mesh,f,aa(1,1),aa(1,2),aa(1,3),aa(1,4),n,key,d0,
     &                  dn)
      a=a+(0.d0,1.0)*aa(:,1)
      b=b+(0.d0,1.0)*aa(:,2)
      c=c+(0.d0,1.0)*aa(:,3)
      d=d+(0.d0,1.0)*aa(:,4)
	deallocate(f,aa)
      end