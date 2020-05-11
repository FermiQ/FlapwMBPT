      subroutine tensor_pack(t6,t)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: t(3,3)
      complex*16, intent(out) :: t6(6)
      integer :: it,i,j
      real*8 :: sq2,dir(3)
      sq2=sqrt(2.d0)/2.d0
      t6=(0.d0,0.d0)
      do it=1,6
        if(it==1) dir=(/1.d0,0.d0,0.d0/)
        if(it==2) dir=(/0.d0,1.d0,0.d0/)
        if(it==3) dir=(/0.d0,0.d0,1.d0/)
        if(it==4) dir=(/sq2,sq2,0.d0/)
        if(it==5) dir=(/sq2,0.d0,sq2/)
        if(it==6) dir=(/0.d0,sq2,sq2/)
        do j=1,3
          do i=1,3
            t6(it)=t6(it)+dir(i)*t(i,j)*dir(j)
          enddo
        enddo
      enddo
      end

      subroutine tensor_pack_r(t6,t)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) :: t(3,3)
      real*8, intent(out) :: t6(6)
      integer :: it,i,j
      real*8 :: sq2,dir(3)
      sq2=sqrt(2.d0)/2.d0
      t6=0.d0
      do it=1,6
        if(it==1) dir=(/1.d0,0.d0,0.d0/)
        if(it==2) dir=(/0.d0,1.d0,0.d0/)
        if(it==3) dir=(/0.d0,0.d0,1.d0/)
        if(it==4) dir=(/sq2,sq2,0.d0/)
        if(it==5) dir=(/sq2,0.d0,sq2/)
        if(it==6) dir=(/0.d0,sq2,sq2/)
        do j=1,3
          do i=1,3
            t6(it)=t6(it)+dir(i)*t(i,j)*dir(j)
          enddo
        enddo
      enddo
      end
