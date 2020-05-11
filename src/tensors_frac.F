      subroutine tensors_frac(t,t1,tf)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: t(3,3),t1(3,3)
      complex*16, intent(out) :: tf(3,3)
      integer :: i
      complex*16 :: t6(6),t16(6)
      call tensor_pack(t6,t)
      call tensor_pack(t16,t1)
      do i=1,6
        t16(i)=t6(i)/t16(i)
      enddo
      call tensor_unpack(t16,tf)
      end


      subroutine tensors_frac_r(t,t1,tf)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) :: t(3,3),t1(3,3)
      real*8, intent(out) :: tf(3,3)
      integer :: i
      real*8 :: t6(6),t16(6)
      call tensor_pack_r(t6,t)
      call tensor_pack_r(t16,t1)
      do i=1,6
        t16(i)=t6(i)/t16(i)
      enddo
      call tensor_unpack_r(t16,tf)
      end
