      subroutine eig_val_gener_zhegvx(n,n0,a,b,new,da)
      implicit none
      integer, intent(in) :: n,n0
      complex*16, intent(inout) :: a(n0,n0),b(n0,n0)
      integer, intent(out) :: new
      real*8, intent(out) :: da(n)
      integer :: info,lwork
      real*8 :: abstol
      integer, allocatable :: iwork(:),ifail(:)
      real*8, allocatable :: rwork(:)
      complex*16, allocatable :: z(:,:),work(:)
      abstol=1.d-50
      lwork=2*max(1,2*n-1)
      allocate(work(lwork),rwork(7*n),iwork(5*n),ifail(n))
      allocate(z(n,n))
      call zhegvx(1,'V','I','U',n,a,n0,b,n0,0.d0,1.d0,1,n,abstol,new,
     &            da,z,n,work,lwork,rwork,iwork,ifail,info)
      call zlacpy(' ',n,n,z,n,a,n0)
      deallocate(z,work,rwork,iwork,ifail)
      end
