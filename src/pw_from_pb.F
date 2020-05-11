      subroutine pw_from_pb(ind_k,ndim,ndim0,pb,pw,key)
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
      complex*16, intent(in) :: pb(ndim0,ndim0)
      complex*16, intent(out) :: pw
      complex*16 :: zdotc
      complex*16, allocatable :: tmp(:)
      if(nplw_gw==0) return
      allocate(tmp(ndim))
      if(key==1) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &           pw_pb(1,ind_k),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_pb(1,ind_k),1,tmp,1)
      else if(key==2) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &           pw_pb_tild(1,ind_k),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_pb(1,ind_k),1,tmp,1)
      else if(key==3) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &             pw_pb_tild(1,ind_k),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_pb_tild(1,ind_k),1,tmp,1)
      endif
      deallocate(tmp)
      end
      
      
      subroutine pw_from_pb_line(ndim,pb,pw,key,pw_tr,pwt_tr)
c     key = 1 - for V, W
c     key = 2 - for E^-1 and E
c     key = 3 - for P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ndim,key
      complex*16, intent(in) :: pb(n_pbtot,n_pbtot),pw_tr(n_pbtot),
     &                          pwt_tr(n_pbtot)
      real*8, intent(out) :: pw
      complex*16 :: zdotc
      complex*16, allocatable :: tmp(:)
      allocate(tmp(ndim))
      if(key==1) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,n_pbtot,pw_tr,1,
     &           (0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_tr,1,tmp,1)
      else if(key==2) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,n_pbtot,pwt_tr,1,
     &             (0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_tr,1,tmp,1)
      else if(key==3) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,n_pbtot,pwt_tr,1,
     &           (0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pwt_tr,1,tmp,1)
      endif
      deallocate(tmp)
      end


      subroutine pw_from_pb_red(iq,ndim,ndim0,pb,pw,key)
c     key = 1 - for V, W
c     key = 2 - for E^-1 and E
c     key = 3 - for P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ndim,ndim0,iq,key
      complex*16, intent(in) :: pb(ndim0,ndim0)
      complex*16, intent(out) :: pw
      complex*16 :: zdotc
      complex*16, allocatable :: tmp(:)
      if(nplw_gw_red==0) return
      allocate(tmp(ndim))
      if(key==1) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &           pw_pb_red(1,iq),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_pb_red(1,iq),1,tmp,1)
      else if(key==2) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &           pw_pb_tild_red(1,iq),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_pb_red(1,iq),1,tmp,1)
      else if(key==3) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &             pw_pb_tild_red(1,iq),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_pb_tild_red(1,iq),1,tmp,1)
      endif
      deallocate(tmp)
      end
      
      
      subroutine pw_from_pb_long(n0,n,pb,pw,ind,key)
c   key = 1 - for V, W
c   key = 2 - for E^-1 and E
c   key = 3 - for P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n0,n,ind,key
      complex*16, intent(in) :: pb(n0,n0)
      real*8, intent(out) :: pw
      complex*16 :: zdotc
      complex*16, allocatable :: tmp(:)
      allocate(tmp(n))
      if(key==1) then
        call zgemv('t',n,n,(1.d0,0.d0),pb,n0,m_pw_long(1,ind),1,
     &           (0.d0,0.d0),tmp,1)
        pw=zdotc(n,m_pw_long(1,ind),1,tmp,1)
      else if(key==2) then
        call zgemv('t',n,n,(1.d0,0.d0),pb,n0,m_pw_tild_long(1,ind),1,
     &           (0.d0,0.d0),tmp,1)
        pw=zdotc(n,m_pw_long(1,ind),1,tmp,1)
      else if(key==3) then
        call zgemv('t',n,n,(1.d0,0.d0),pb,n0,m_pw_tild_long(1,ind),1,
     &           (0.d0,0.d0),tmp,1)
        pw=zdotc(n,m_pw_tild_long(1,ind),1,tmp,1)
      endif
      deallocate(tmp)
      end


      subroutine pw_from_pb_1q(ndim,pb,pw,key,pw_tr,pwt_tr)
c     key = 1 - for V, W
c     key = 2 - for E^-1 and E
c     key = 3 - for P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ndim,key
      complex*16, intent(in) :: pb(ndim,ndim),pw_tr(ndim),pwt_tr(ndim)
      complex*16, intent(out) :: pw
      complex*16 :: zdotc
      complex*16, allocatable :: tmp(:)
      allocate(tmp(ndim))
      if(key==1) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim,pw_tr,1,
     &           (0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_tr,1,tmp,1)
      else if(key==2) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim,pwt_tr,1,
     &             (0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_tr,1,tmp,1)
      else if(key==3) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim,pwt_tr,1,
     &           (0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pwt_tr,1,tmp,1)
      endif
      deallocate(tmp)
      end
