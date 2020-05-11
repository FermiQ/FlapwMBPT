      subroutine spline_smoothing(n,mu,hh,rrr,ff,sg,d,l1,l2,y,as)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: mu,rrr(n),ff(n),sg(0:n),d(n-1),l1(n-2),
     &                      l2(n-3),y(0:n),hh(n)
      real*8, intent(out) :: as(n+1,4)
      integer :: i
      real*8, allocatable :: zz(:)
      allocate(zz(n-1))
      do i=1,n-1
        zz(i)=rrr(i)*y(i-1)+ff(i)*y(i)+rrr(i+1)*y(i+1)
      enddo
      as(1,3)=0.d0
      as(n+1,3)=0.d0
      call pentadiag_solve(n-1,d,l1,l2,zz,as(2,3))
      deallocate(zz)
c ------- Reconstruct A(1:N+1) ------------------------------
      as(1,1)=rrr(1)*as(2,3)
      as(2,1)=ff(1)*as(2,3)+rrr(2)*as(3,3)
      do i=3,n-1
        as(i,1)=rrr(i-1)*as(i-1,3)+ff(i-1)*as(i,3)+rrr(i)*as(i+1,3)
      enddo
      as(n,1)=rrr(n-1)*as(n-1,3)+ff(n-1)*as(n,3)
      as(n+1,1)=rrr(n)*as(n,3)
      do i=1,n+1
        as(i,1)=y(i-1)-mu*sg(i-1)*as(i,1)
      enddo
c ------- Reconstruct B(1:N) ------------------------------
      do i=1,n
        as(i,2)=(as(i+1,1)-as(i,1))/hh(i)
     &         -(as(i+1,3)+2.d0*as(i,3))*hh(i)/3.d0
      enddo
c ------- Reconstruct D(1:N) ------------------------------
      do i=1,n
        as(i,4)=(as(i+1,3)-as(i,3))/hh(i)/3.d0
      enddo
      end