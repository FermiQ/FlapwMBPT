      subroutine tensor_inverse(t,ti)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: t(3,3)
      complex*16, intent(out) :: ti(3,3)
      integer :: i
      complex*16 :: t6(6)
      call tensor_pack(t6,t)
      do i=1,6
        t6(i)=(1.d0,0.d0)/t6(i)
      enddo
      call tensor_unpack(t6,ti)
      end


      subroutine tensor_inverse_r(t,ti)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) :: t(3,3)
      real*8, intent(out) :: ti(3,3)
      integer :: i
      real*8 :: t6(6)
      call tensor_pack_r(t6,t)
      do i=1,6
        t6(i)=1.d0/t6(i)
      enddo
      call tensor_unpack_r(t6,ti)
      end
