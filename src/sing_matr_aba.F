      subroutine sing_matr_aba(n,b00,b01,b10,b11,a00,a01,a11,c00,c01,
     &                         c11)
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
      complex*16, intent(in) :: a00(3,3),a01(3,n),a11(n,n),b00(3,3),
     &                          b01(3,n),b10(3,n),b11(n,n)
      complex*16, intent(out) :: c00(3,3),c01(3,n),c11(n,n)
      integer :: i,ii,j,jj
      complex*16, allocatable :: d00(:,:),d01(:,:),d10(:,:),d11(:,:)
      allocate(d00(3,3))
      allocate(d01(3,n))
      allocate(d10(n,3))
      allocate(d11(n,n))
      d00=matmul(b00,a00)
      do jj=1,3
        do ii=1,3
          do i=1,n
            d00(ii,jj)=d00(ii,jj)+b01(ii,i)*conjg(a01(jj,i))
          enddo
        enddo
      enddo
      d01=(0.d0,0.d0)
      do j=1,n
        do ii=1,3
          do jj=1,3
            d01(ii,j)=d01(ii,j)+b00(ii,jj)*a01(jj,j)
          enddo
          do i=1,n
            d01(ii,j)=d01(ii,j)+b01(ii,i)*a11(i,j)
          enddo
        enddo
      enddo
      d10=(0.d0,0.d0)
      do jj=1,3
        do i=1,n
          do ii=1,3
            d10(i,jj)=d10(i,jj)+b10(ii,i)*a00(ii,jj)
          enddo
          do j=1,n
            d10(i,jj)=d10(i,jj)+b11(i,j)*conjg(a01(jj,j))
          enddo
        enddo
      enddo
      d11=matmul(b11,a11)
      do j=1,n
        do i=1,n
          do ii=1,3
            d11(i,j)=d11(i,j)+b10(ii,i)*a01(ii,j)
          enddo
        enddo
      enddo
c ------------------------------------------------------------------
      c00=matmul(d00,b00)
      do jj=1,3
        do ii=1,3
          do i=1,n
            c00(ii,jj)=c00(ii,jj)+d01(ii,i)*conjg(b01(jj,i))
          enddo
        enddo
      enddo
      c01=(0.d0,0.d0)
      do j=1,n
        do ii=1,3
          do jj=1,3
            c01(ii,j)=c01(ii,j)+d00(ii,jj)*b01(jj,j)
          enddo
          do i=1,n
            c01(ii,j)=c01(ii,j)+d01(ii,i)*b11(i,j)
          enddo
        enddo
      enddo
      c11=matmul(d11,b11)
      do j=1,n
        do i=1,n
          do ii=1,3
            c11(i,j)=c11(i,j)+d10(ii,i)*b01(ii,j)
          enddo
        enddo
      enddo
      deallocate(d00,d01,d10,d11)
      end
