      subroutine pentadiag_fact(n,a0,a1,a2,d,l1,l2)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: a0(n),a1(n-1),a2(n-2)
      real*8, intent(out) :: d(n),l1(n-1),l2(n-2)
      integer :: i
      d(1)=a0(1)
      l1(1)=a1(1)/d(1)
      l2(1)=a2(1)/d(1)
      d(2)=a0(2)-l1(1)**2*d(1)
      do i=2,n-2
        l1(i)=(a1(i)-l1(i-1)*l2(i-1)*d(i-1))/d(i)
        l2(i)=a2(i)/d(i)
        d(i+1)=a0(i+1)-l2(i-1)**2*d(i-1)-l1(i)**2*d(i)
      enddo
      l1(n-1)=(a1(n-1)-l1(n-2)*l2(n-2)*d(n-2))/d(n-1)
      d(n)=a0(n)-l2(n-2)**2*d(n-2)-l1(n-1)**2*d(n-1)
      end