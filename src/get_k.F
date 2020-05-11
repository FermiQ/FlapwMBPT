      subroutine get_k(ierror)
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
	  integer :: ierror
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call sbzmesh(pnt,wgt,npnt,ierror)
      if(maswrk) write(iun,1030)npnt
1030  format(1x,i4,' k-points generated')
      end
