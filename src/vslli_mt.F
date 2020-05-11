      subroutine vslli_mt(key,key1,vhmt,vxcmt)
c	  key=0 - only Hartree part
c	  key=1 - Only XC part
c	  key=2 - Hartree + XC
c	  key1=1 - L=0 component is taken into account
c	  key1=2 - L=0 component is not taken into account
	  use atom_mod
	  use etot_mod
	  use manager_mod
	  use solid_mod
      implicit none
	  integer, intent(in) :: key,key1
	  real*8, intent(in) :: vhmt(maxmt/nspin_0),vxcmt(maxmt)
	  integer :: ispin,mt
	  real*8, allocatable :: v(:)
c ---------------------------------------------------------------------
      allocate(v(maxmt/nspin_0))
	  if(key==0) then
	    v=vhmt
	  endif
	  do ispin=1,nspin
	    mt=indmt(1,1,ispin)
	    if(key==1) then
	      v=vxcmt(mt:mt+maxmt/nspin_0-1)
	    else if(key==2) then
	      v=vhmt+vxcmt(mt:mt+maxmt/nspin_0-1)
	    endif
	    call v_me(ispin,key1,key,v,0)
	  enddo
c ----------------------------------------------------------------	  
      call sym_me
	  deallocate(v)
      call timel('******* vslli_mt finished **********')
      end
