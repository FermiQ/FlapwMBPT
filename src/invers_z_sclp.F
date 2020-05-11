#ifdef SCLP
      subroutine invers_z_sclp(n,nrow,ncol,nblock,desc,a)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol,desc(9),nblock
      complex*16, intent(inout) :: a(nrow,ncol)
      integer :: info,lwork,liwork
      integer, allocatable :: ipiv(:),iwork(:)
      complex*16,allocatable :: wrk(:)
      allocate(ipiv(nrow+nblock))
      lwork=(nrow+nblock)*nblock
      allocate(wrk(lwork))
      liwork=lwork
      allocate(iwork(liwork))
      call pzgetrf(n,n,a,1,1,desc,ipiv,info)
      call pzgetri(n,a,1,1,desc,ipiv,wrk,lwork,iwork,liwork,info)
      deallocate(ipiv,wrk,iwork)
      end


      subroutine invers_r_sclp(n,nrow,ncol,nblock,desc,a)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol,desc(9),nblock
      real*8, intent(inout) :: a(nrow,ncol)
      integer :: info,lwork,liwork
      integer, allocatable :: ipiv(:),iwork(:)
      real*8,allocatable :: wrk(:)
      allocate(ipiv(nrow+nblock))
      lwork=(nrow+nblock)*nblock
      allocate(wrk(lwork))
      liwork=lwork
      allocate(iwork(liwork))
      call pdgetrf(n,n,a,1,1,desc,ipiv,info)
      call pdgetri(n,a,1,1,desc,ipiv,wrk,lwork,iwork,liwork,info)
      deallocate(ipiv,wrk,iwork)
      end
#endif
