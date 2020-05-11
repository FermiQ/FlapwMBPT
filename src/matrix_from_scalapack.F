      subroutine matrix_from_scalapack(n0,a,nrow,ncol,aloc,refrow,
     &                                 refcol)
      use parallel_mod
      implicit none
      integer, intent(in) :: n0,nrow,ncol,refrow(nrow),refcol(ncol)
      complex*16, intent(in) :: aloc(nrow,ncol)
      complex*16, intent(out) :: a(n0,n0)
      integer :: i,j,ii,jj
      a=(0.d0,0.d0)
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          a(ii,jj)=aloc(i,j)
        enddo
      enddo
      call dgop(a,2*n0*n0,'  +',comm_b)
      end

      subroutine matrix_from_scalapack_r(n0,a,nrow,ncol,aloc,refrow,
     &                                   refcol)
      use parallel_mod
      implicit none
      integer, intent(in) :: n0,nrow,ncol,refrow(nrow),refcol(ncol)
      real*8, intent(in) :: aloc(nrow,ncol)
      real*8, intent(out) :: a(n0,n0)
      integer :: i,j,ii,jj
      a=0.d0
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          a(ii,jj)=aloc(i,j)
        enddo
      enddo
      call dgop(a,n0*n0,'  +',comm_b)
      end
