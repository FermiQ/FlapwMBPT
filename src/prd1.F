      subroutine prd1(a,b,c)
	  integer :: i2,i3
      real*8 :: a(3),b(3,3),c(3)
      do i2=1,3
        c(i2)=0.d0
        do i3=1,3
          c(i2)=c(i2)+a(i3)*b(i3,i2)
        enddo
      enddo
      end
