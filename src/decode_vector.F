      subroutine decode_vector(t,i1,i2,r3)
      character*80, intent(in) :: t
      integer, intent(in) :: i1,i2
      real*8, intent(out) :: r3(3)
      integer :: j1,j2
      j1=i1
      j2=j1+index(t(i1:i2),',')-2
      call decode_number(t,j1,j2,r3(1))
      j1=j2+2
      j2=j1+index(t(j1:i2),',')-2
      call decode_number(t,j1,j2,r3(2))
      j1=j2+2
      call decode_number(t,j1,i2,r3(3))
      end
