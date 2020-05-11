      SUBROUTINE ro_int_g(correlation,rhoi)
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      IMPLICIT none
      logical, intent(in) :: correlation
      complex*16, intent(out) :: rhoi(nplwro,nspin)
	  integer :: ispin,k,ind_k
	  logical :: jreal,jimag
	  rhoi=(0.d0,0.d0)
      DO ISPIN=1,NSPIN
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
	      call sumi_gw(correlation,ind_k,k,ispin,rhoi)
        enddo 
        if(nproc_k/=1) then
	      call DGOP(rhoi(1,ispin),2*nplwro,'  +',comm_k)
          if(magn.eq.2) then
	        call DGOP(spintr,6*nplwro,'  +',comm_k)
	      endif
	    endif
c ------- Symmetrization ----------------------------------------------
        call symscal(rhoi(1,ispin))
      enddo   !!! over ispin
      if(magn.eq.2) then
        jreal=.false.
        jimag=.false.
        if(inv.lt.2) jreal=.true.
        if(inv.ne.1) jimag=.true.
        call symvec(spintr,jreal,jimag)
      endif  !!! for magn = 2 only
      END
