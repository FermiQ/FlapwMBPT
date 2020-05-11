      subroutine spline_interp_hmg(f1,f2,n1,n2,h1,h2)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n1,n2
      real*8, intent(in) :: h1,h2,f1(0:n1)
      real*8, intent(out) :: f2(0:n2)
      integer :: i,i0
      real*8 :: x,dx
      real*8, allocatable :: a(:),b(:),c(:),d(:)
      allocate(a(n1))
      allocate(b(n1))
      allocate(c(n1))
      allocate(d(n1))
      call spline_hmg(h1,f1,a,b,c,d,n1,0,0.d0,0.d0)
      do i=0,n2
        x=h2*i
        i0=x/h1+1
        i0=max(i0,0)
        i0=min(i0,n1)
        f2(i)=a(i0)+b(i)*dx+c(i0)*dx*dx+d(i0)*dx*dx*dx
      enddo
      deallocate(a,b,c,d)
      end