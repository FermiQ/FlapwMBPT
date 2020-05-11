      subroutine theta_xc_pw(key)
c     key = 0 - result is provided in bosonic 4-vector
c     key = 1 - the result is provided in spin-spin 2*2 form
      use atom_mod
      use manager_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      allocate(theta_mt(maxmt/nspin,n_ixc_0))
      call theta_mt_xc(key)
      if(allocated(theta_int)) deallocate(theta_int)
      allocate(theta_int(nplwro,n_ixc_0))
      call theta_int_xc(key)
      call timel('******* THETA_XC_PW finished *******')
      end
