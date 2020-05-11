      subroutine c_bessel(lam,r,lmax,cl)
      implicit none
      integer, intent(in) :: lmax
      real*8, intent(in) :: r,lam
      real*8, intent(out) :: cl(0:lmax)
      integer :: l
      real*8 :: x
      real*8, allocatable :: gam(:)
      allocate(gam(0:lmax))
      gam(0)=1.d0
      gam(1)=1.d0+lam*r
      x=(lam*r)**2
      do l=1,lmax-1
        gam(l+1)=(l+l+1.d0)*gam(l)+x*gam(l-1)
      enddo
      x=exp(-lam*r)
      cl=gam
      if(r>1.d-12) then
        do l=0,lmax
          cl(l)=x*gam(l)/r**(l+1)
        enddo
      endif
      deallocate(gam)
      end
