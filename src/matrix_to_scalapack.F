      subroutine matrix_to_scalapack(n0,a,nrow,ncol,aloc,refrow,refcol)
      use parallel_mod
      implicit none
      integer, intent(in) :: n0,nrow,ncol,refrow(nrow),refcol(ncol)
      complex*16, intent(in) :: a(n0,n0)
      complex*16, intent(out) :: aloc(nrow,ncol)
      integer :: i,j,ii,jj
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          aloc(i,j)=a(ii,jj)
        enddo
      enddo
      end

      subroutine matrix_to_scalapack_r(n0,a,nrow,ncol,aloc,refrow,
     &                                 refcol)
      use parallel_mod
      implicit none
      integer, intent(in) :: n0,nrow,ncol,refrow(nrow),refcol(ncol)
      real*8, intent(in) :: a(n0,n0)
      real*8, intent(out) :: aloc(nrow,ncol)
      integer :: i,j,ii,jj
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          aloc(i,j)=a(ii,jj)
        enddo
      enddo
      end
