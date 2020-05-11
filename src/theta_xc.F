      subroutine theta_xc(key)
c     key = 0 - result is provided in bosonic 4-vector
c     key = 1 - the result is provided in spin-spin 2*2 form
	  use atom_mod
	  use manager_mod
	  use solid_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: key
      integer :: i
      real*8, allocatable :: tmp(:,:,:)
	  allocate(theta_mt(maxmt/nspin,n_ixc_0))
      call theta_mt_xc(key)
      if(allocated(thet_me)) deallocate(thet_me)
      allocate(thet_me(n_pbmtm,n_pbmtm,natom,n_ixc_0))
      call theta_mtm
      deallocate(theta_mt)
      if(ncmpl==1) then
        allocate(tmp(n_pbmtm,n_pbmtm,natom))
        do i=1,n_ixc_0
          call theta_r_from_c(thet_me(1,1,1,i),tmp)
          thet_me(:,:,:,i)=tmp
        enddo
        deallocate(tmp)
      endif
      call timel('******* THETA_XC_MT finished *******')
      if(allocated(theta_int)) deallocate(theta_int)
	  allocate(theta_int(nplwdif_pb,n_ixc_0))
      call theta_int_xc(key)
      call timel('******* THETA_XC finished **********')
      end
