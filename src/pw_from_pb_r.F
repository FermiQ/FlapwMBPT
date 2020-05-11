      subroutine pw_from_pb_r(ind_k,ndim,ndim0,pb,pw,key)
c     key = 1 - for V, W
c     key = 2 - for E^-1 and E
c     key = 3 - for P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ndim,ndim0,ind_k,key
      real*8, intent(in) :: pb(ndim0,ndim0)
      real*8, intent(out) :: pw
      real*8, allocatable :: tmp(:)
      if(nplw_gw==0) return
      allocate(tmp(ndim))
      if(key==1) then
        call dgemv('n',ndim,ndim,1.d0,pb,ndim0,pw_pb_r(1,ind_k),1,0.d0,
     &             tmp,1)
        pw=dot_product(pw_pb_r(1:ndim,ind_k),tmp)
      else if(key==2) then
        call dgemv('n',ndim,ndim,1.d0,pb,ndim0,pw_pb_tild_r(1,ind_k),1,
     &             0.d0,tmp,1)
        pw=dot_product(pw_pb_r(1:ndim,ind_k),tmp)
      else if(key==3) then
        call dgemv('n',ndim,ndim,1.d0,pb,ndim0,pw_pb_tild_r(1,ind_k),1,
     &             0.d0,tmp,1)
        pw=dot_product(pw_pb_tild_r(1:ndim,ind_k),tmp)
      endif
      deallocate(tmp)
      end
