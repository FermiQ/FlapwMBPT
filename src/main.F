      Program Main
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
c***********************************************************
!              The main module                             *
!***********************************************************
      implicit none
      call beging
      mem_using=0.d0
      mem_max=0.d0
      mem_where='qft'
      call init
      call set_control
      if(maswrk) write(iun,*)' Number of processors is ',nproc
      call start_par_mb
! --------------------------------------------------------------
      obietto='HEG'
      if(rs_el_gas<0.d0) then
        obietto='SOL'
        rs_el_gas=-rs_el_gas
      endif
      if(obietto=='SOL') then
c        call v0_b0_finder
        call tau_frequency_meshes_spl
        call solid
      else if(obietto=='HEG') then
        iheg=.true.
        call heg
      endif
!     ***********************************************************
      call ending
      end
