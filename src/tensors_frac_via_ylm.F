      subroutine tensors_frac_via_ylm(t,t1,tf,lmax)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmax
      complex*16, intent(in) :: t(3,3),t1(3,3)
      complex*16, intent(out) :: tf(3,3)
      integer :: lmmax
      complex*16 :: ty(9),t1y(9)
      complex*16, allocatable :: ry(:)
      lmmax=(lmax+1)**2
      call tensor_y_exp(t,ty)
      call tensor_y_exp(t1,t1y)
      allocate(ry(lmmax))
      ry=(0.d0,0.d0)
      ry(1:9)=ty
      call ratio_ylm(lmmax,t1y,ry)
      call tensor_from_ylm(lmax,tf,ry)
      deallocate(ry)
      end
