      SUBROUTINE ROfull(correlation,rhom,rhoi)
**********************************************************
* The calculation of the  valence density                *
**********************************************************
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      IMPLICIT none
      logical, intent(in) :: correlation
      real*8, intent(out) :: rhom(maxmt)
      complex*16, intent(out) :: rhoi(nplwro,nspin)
      logical :: jreal,jimag
	  integer :: keym
	  complex*16, allocatable :: gloc(:,:,:,:)
      if(magn.eq.2) spintr=(0.d0,0.d0)
	  call ro_int_g(correlation,rhoi)
      call timel('****** Rho_Int finished ************')
c ------------------------------------------------------------------
	  rhom=0.d0
      if(magn.eq.2) spmt=0.d0
      allocate(gloc(maxel,maxel,natom,nspin))
      call g_mt_from_g_full(correlation,gloc)
      call get_t(gloc)
      deallocate(gloc)
	  call ro_val(rhom)
	  call renr_val(rhom,rhoi)
      call timel('****** Rho_Mt finished *************')
	  call add_cor(1)
      if(itest>0) then
        keym=0
        if(nspin.eq.2) keym=2
	    text=' Matching : New Density!  '
        call matching(keym,rhom,rhoi)
        if(nspin.eq.2) then
	      text=' Matching : New Density of Spin Moment - Z!'
          call matching(3,rhom,rhoi)
	    endif
	    if(magn.eq.2) then
          jreal=.false.
          jimag=.false.
          if(inv.lt.2) jreal=.true.
          if(inv.ne.1) jimag=.true.
	      text=' Matching : New Density of Spin Moment!'
          call matchav(spmt,spintr,jreal,jimag)
        endif
      endif  !!! over itest
	  call mag_mom
      call timel('****** Rofull finished *************')
      END
