#ifdef SCLP
      subroutine eig_val_only_sclp(n,a,da,nrow,ncol,desc)
      implicit none
      integer, intent(in) :: n,nrow,ncol,desc(9)
      complex*16, intent(inout) :: a(nrow,ncol)
      real*8, intent(out) :: da(n)
      integer :: lwork,info,new,nz,id,lrwork,liwork,nn
      real*8 :: vl,vu,abstol
      complex*16 :: z
      integer, allocatable :: iwork(:),ifail(:)
      real*8, allocatable :: rwork(:)
      complex*16, allocatable :: work(:)
      abstol=1.d-70
      vl=-1.d10
      vu=1.d10
      lwork=n+max(3*64,64*(max(nrow,ncol)+1))
      lwork=2*lwork
      nn=max(n,64)
      lrwork=5*nn+4*n
      nn=max(n,nrow*ncol+1,4)
      liwork=6*nn
      allocate(work(lwork),rwork(lrwork),iwork(liwork),ifail(n))
      call pzheevx('N','A','L',n,a,1,1,desc,vl,vu,1,1,abstol,new,nz,da,
     &             1.d0,z,1,1,desc,work,lwork,rwork,lrwork,iwork,liwork,
     &             ifail,id,z,info)
      deallocate(work,rwork,iwork,ifail)
      end

      subroutine eig_val_only_sclp_r(n,a,da,nrow,ncol,desc)
      implicit none
      integer, intent(in) :: n,nrow,ncol,desc(9)
      real*8, intent(inout) :: a(nrow,ncol)
      real*8, intent(out) :: da(n)
      integer :: lwork,info,new,nz,id,liwork,nn
      real*8 :: vl,vu,abstol
      complex*16 :: z
      integer, allocatable :: iwork(:),ifail(:)
      real*8, allocatable :: work(:)
      abstol=1.d-70
      vl=-1.d10
      vu=1.d10
      lwork=n+max(3*64,64*(max(nrow,ncol)+1))
      lwork=2*lwork
      nn=max(n,nrow*ncol+1,4)
      liwork=6*nn
      allocate(work(lwork),iwork(liwork),ifail(n))
      call pdsyevx('N','A','L',n,a,1,1,desc,vl,vu,1,1,abstol,new,nz,da,
     &             1.d0,z,1,1,desc,work,lwork,iwork,liwork,ifail,id,z,
     &             info)
      deallocate(work,iwork,ifail)
      end
#endif
