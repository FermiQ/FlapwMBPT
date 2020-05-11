      subroutine pw_from_pb_sclp(ind_k,nrow,ncol,refrow,refcol,pb,pw,
     &                           key)
c     key = 1 - for V, W
c     key = 2 - for E^-1 and E
c     key = 3 - for P, R
      use atom_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: nrow,ncol,ind_k,key,refrow(nrow),
     &                       refcol(ncol)
      complex*16, intent(in) :: pb(nrow,ncol)
      complex*16, intent(out) :: pw
      integer :: i,j,ii,jj
      complex*16 :: cc
      pw=(0.d0,0.d0)
      if(key==1) then
        do j=1,ncol
          jj=refcol(j)
          cc=(0.d0,0.d0)
          do i=1,nrow
            ii=refrow(i)
            cc=cc+conjg(pw_pb(ii,ind_k))*pb(i,j)
          enddo
          pw=pw+cc*pw_pb(jj,ind_k)
        enddo
      else if(key==2) then
        do j=1,ncol
          jj=refcol(j)
          cc=(0.d0,0.d0)
          do i=1,nrow
            ii=refrow(i)
            cc=cc+conjg(pw_pb(ii,ind_k))*pb(i,j)
          enddo
          pw=pw+cc*pw_pb_tild(jj,ind_k)
        enddo
      else if(key==3) then
        do j=1,ncol
          jj=refcol(j)
          cc=(0.d0,0.d0)
          do i=1,nrow
            ii=refrow(i)
            cc=cc+conjg(pw_pb_tild(ii,ind_k))*pb(i,j)
          enddo
          pw=pw+cc*pw_pb_tild(jj,ind_k)
        enddo
      endif
      call dgop(pw,2,'  +',comm_b)
      end
