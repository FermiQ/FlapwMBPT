      subroutine size_shift_par(n,npr,ndim_pr,n_mpi_pr)
      implicit none
	  integer, intent(in) :: n,npr
	  integer, intent(out) :: ndim_pr(npr),n_mpi_pr(npr)
	  integer :: ish,irest,ind,i
	  ish=n/npr
	  irest=n-ish*npr
	  ndim_pr(1:irest)=ish+1
	  ndim_pr(irest+1:npr)=ish
	  ind=0
	  do i=1,npr
	    n_mpi_pr(i)=ind
	    ind=ind+ndim_pr(i)
	  enddo
	  end
