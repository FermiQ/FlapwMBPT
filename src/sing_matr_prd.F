      subroutine sing_matr_prd(a00,a01,a10,a11,b00,b01,b10,b11,c00,c01,
     &                         c10,c11,n)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      complex*16, intent(in) :: a00(3,3),a01(3,n),a10(3,n),a11(n,n),
     &                          b00(3,3),b01(3,n),b10(3,n),b11(n,n)
      complex*16, intent(out) :: c00(3,3),c01(3,n),c10(3,n),c11(n,n)
      integer :: i,ii,j,jj
      c00=matmul(a00,b00)
      do jj=1,3
        do ii=1,3
          do i=1,n
            c00(ii,jj)=c00(ii,jj)+a01(ii,i)*b10(jj,i)
          enddo
        enddo
      enddo
      c01=(0.d0,0.d0)
      do j=1,n
        do ii=1,3
          do jj=1,3
            c01(ii,j)=c01(ii,j)+a00(ii,jj)*b01(jj,j)
          enddo
          do i=1,n
            c01(ii,j)=c01(ii,j)+a01(ii,i)*b11(i,j)
          enddo
        enddo
      enddo
      c10=(0.d0,0.d0)
      do jj=1,3
        do i=1,n
          do ii=1,3
            c10(jj,i)=c10(jj,i)+a10(ii,i)*b00(ii,jj)
          enddo
          do j=1,n
            c10(jj,i)=c10(jj,i)+a11(i,j)*b10(jj,j)
          enddo
        enddo
      enddo
      c11=matmul(a11,b11)
      do j=1,n
        do i=1,n
          do ii=1,3
            c11(i,j)=c11(i,j)+a10(ii,i)*b01(ii,j)
          enddo
        enddo
      enddo
      end
