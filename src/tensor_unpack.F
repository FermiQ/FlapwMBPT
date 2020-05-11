      subroutine tensor_unpack(t6,t)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: t6(6)
      complex*16, intent(out) :: t(3,3)
      t(1,1)=t6(1)
      t(2,2)=t6(2)
      t(3,3)=t6(3)
      t(1,2)=t6(4)-(t(1,1)+t(2,2))/2
      t(1,3)=t6(5)-(t(1,1)+t(3,3))/2
      t(2,3)=t6(6)-(t(2,2)+t(3,3))/2
      t(2,1)=t(1,2)
      t(3,1)=t(1,3)
      t(3,2)=t(2,3)
      end


      subroutine tensor_unpack_r(t6,t)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) :: t6(6)
      real*8, intent(out) :: t(3,3)
      t(1,1)=t6(1)
      t(2,2)=t6(2)
      t(3,3)=t6(3)
      t(1,2)=t6(4)-(t(1,1)+t(2,2))/2
      t(1,3)=t6(5)-(t(1,1)+t(3,3))/2
      t(2,3)=t6(6)-(t(2,2)+t(3,3))/2
      t(2,1)=t(1,2)
      t(3,1)=t(1,3)
      t(3,2)=t(2,3)
      end
