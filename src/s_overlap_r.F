      subroutine s_overlap_r(key,s,n0,n,iref)
c     Calculates interstitial overlap matrix for plane waves ------
c     IREF - reference of the plane waves to GBS array ------------
c     key = 0  Just overlap matrix
c     key = 1  Its inverse
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: key,n,n0,iref(n0)
      real*8, intent(out) :: s(n0,n0)
      integer :: j,j0,i,i0,ind,info,i1,i2,i3
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: work(:)
      do j=1,n
        j0=iref(j)
        do i=1,n
          i0=iref(i)
          i1=igbs(1,j0)-igbs(1,i0)
          i2=igbs(2,j0)-igbs(2,i0)
          i3=igbs(3,j0)-igbs(3,i0)
          ind=indplw(i1,i2,i3)
          s(i,j)=sovr(ind)
        enddo
      enddo
      if(key==1) then
        allocate(ipiv(n),work(3*n))
        call dgetrf(n,n,s,n0,ipiv,info)
        call dgetri(n,s,n0,ipiv,work,3*n,info)
        deallocate(ipiv,work)
      endif
      end
