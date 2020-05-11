      subroutine b_bessel(lam,r,lmax,bl)
      implicit none
      integer, intent(in) :: lmax
      real*8, intent(in) :: r,lam
      real*8, intent(out) :: bl(0:lmax)
      integer :: l
      complex*16 :: k,x
      complex*16, allocatable :: bc(:),bcd(:)
      allocate(bc(0:lmax))
      allocate(bcd(0:lmax))
      k=(0.d0,1.d0)*lam
      x=k*r
      call bessc(x,bc,bcd,lmax)
      do l=0,lmax
        bl(l)=bc(l)/k**l
      enddo
      deallocate(bc,bcd)
      end
