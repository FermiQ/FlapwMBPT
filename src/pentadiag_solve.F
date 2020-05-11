      subroutine pentadiag_solve(n,d,l1,l2,zz,c)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: d(n),l1(n-1),l2(n-2),zz(n)
      real*8, intent(out) :: c(n)
      integer :: i
      real*8, allocatable :: x(:)
      allocate(x(n))
      x(1)=zz(1)
      x(2)=zz(2)-l1(1)*x(1)
      do i=3,n
        x(i)=zz(i)-l2(i-2)*x(i-2)-l1(i-1)*x(i-1)
      enddo
      do i=1,n
        x(i)=x(i)/d(i)
      enddo
      c(n)=x(n)
      c(n-1)=x(n-1)-l1(n-1)*c(n)
      do i=n-2,1,-1
        c(i)=x(i)-l1(i)*c(i+1)-l2(i)*c(i+2)
      enddo
      deallocate(x)
      end