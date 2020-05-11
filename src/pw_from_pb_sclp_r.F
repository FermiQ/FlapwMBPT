      subroutine pw_from_pb_sclp_r(ind_k,nrow,ncol,refrow,refcol,pb,pw,
     &                             key)
c     key = 1 - for V, W
c     key = 2 - for E^-1 and E
c     key = 3 - for P, R
      use atom_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: nrow,ncol,ind_k,key,refrow(nrow),
     &                       refcol(ncol)
      real*8, intent(in) :: pb(nrow,ncol)
      real*8, intent(out) :: pw
      integer :: i,j,ii,jj
      real*8 :: cc
      pw=0.d0
      if(key==1) then
        do j=1,ncol
          jj=refcol(j)
          cc=0.d0
          do i=1,nrow
            ii=refrow(i)
            cc=cc+pw_pb_r(ii,ind_k)*pb(i,j)
          enddo
          pw=pw+cc*pw_pb_r(jj,ind_k)
        enddo
      else if(key==2) then
        do j=1,ncol
          jj=refcol(j)
          cc=0.d0
          do i=1,nrow
            ii=refrow(i)
            cc=cc+pw_pb_r(ii,ind_k)*pb(i,j)
          enddo
          pw=pw+cc*pw_pb_tild_r(jj,ind_k)
        enddo
      else if(key==3) then
        do j=1,ncol
          jj=refcol(j)
          cc=0.d0
          do i=1,nrow
            ii=refrow(i)
            cc=cc+pw_pb_tild_r(ii,ind_k)*pb(i,j)
          enddo
          pw=pw+cc*pw_pb_tild_r(jj,ind_k)
        enddo
      endif
      call dgop(pw,1,'  +',comm_b)
      end
