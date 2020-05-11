      subroutine tensor_inverse_via_ylm(t,ti,lmax)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmax
      complex*16, intent(in) :: t(3,3)
      complex*16, intent(out) :: ti(3,3)
      integer :: lmmax
      complex*16 :: ty(9)
      complex*16, allocatable :: ry(:)
      lmmax=(lmax+1)**2
      call tensor_y_exp(t,ty)
      allocate(ry(lmmax))
      ry=(0.d0,0.d0)
      ry(1)=sqrt(4.d0*pi)
      call ratio_ylm(lmmax,ty,ry)
      call tensor_from_ylm(lmax,ti,ry)
      deallocate(ry)
      end
