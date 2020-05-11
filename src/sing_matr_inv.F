      subroutine sing_matr_inv(a00,a01,a11,n)
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
      complex*16, intent(inout) :: a00(3,3),a01(3,n),a11(n,n)
      integer :: i,ii,j,jj
      complex*16 :: a00inv(3,3),a(3,3)
      complex*16, allocatable :: b(:,:),c(:,:)
c ---------- Evaluate a00^(-1) -------------------------------------
      call tensor_inverse(a00,a00inv)
c ---------- Evaluate a11^(-1) -------------------------------------
      allocate(b(n,n))
      b=a11
      call invers_h(n,b,n)
c -------------------------- Final a00 -----------------------------
      a=a00
      do jj=1,3
        do ii=1,3
          do j=1,n
            do i=1,n
c              a(ii,jj)=a(ii,jj)-conjg(a01(ii,i))*b(i,j)*a01(jj,j)
            enddo
          enddo
        enddo
      enddo
      call tensor_inverse(a,a00)
c ------------------------------------------------------------------
      do j=1,n
        do i=1,n
          do jj=1,3
            do ii=1,3
c              a11(i,j)=a11(i,j)-conjg(a01(ii,i))*a00inv(ii,jj)*a01(jj,j)
            enddo
          enddo
        enddo
      enddo
      call invers_h(n,a11,n)
c ------------------------------------------------------------------
      allocate(c(3,n))
      c=(0.d0,0.d0)
      do j=1,n
        do ii=1,3
          do i=1,n
            c(ii,j)=c(ii,j)+a01(ii,i)*a11(i,j)
          enddo
        enddo
      enddo
      a01=(0.d0,0.d0)
      do j=1,n
        do ii=1,3
          do jj=1,3
            a01(ii,j)=a01(ii,j)-a00inv(ii,jj)*c(jj,j)
          enddo
        enddo
      enddo
      deallocate(b,c)
      end
