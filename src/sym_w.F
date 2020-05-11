      subroutine sym_w(iq,iq0,w_0,n_pbt)
	  use atom_mod
	  use solid_mod
      implicit none
      integer, intent(in) :: iq,iq0,n_pbt
	  complex*16, intent(inout) :: w_0(n_pbt,n_pbt)
	  integer :: j,i,ig,npw
	  complex*16, allocatable :: tmp(:,:)
c -------- Copy+Symmetrization of W into another k-point ---------------
	  ig=k_group(iq)
	  npw=nplwgw(iq0)
      if(iq/=iq0) then
	    call sym_w_mt(ig,pnt(1,iq),w_0,n_pbt)
	    allocate(tmp(n_pbmt,nplw_gw))
	    tmp(1:n_pbmt,1:npw)=w_0(1:n_pbmt,n_pbmt+1:n_pbt)
	    call sym_w_mi(iq,ig,tmp)
	    w_0(1:n_pbmt,n_pbmt+1:n_pbt)=tmp(1:n_pbmt,1:npw)
	    deallocate(tmp)
	    allocate(tmp(npw,npw))
	    tmp=w_0(n_pbmt+1:n_pbt,n_pbmt+1:n_pbt)
	    call sym_w_ii(iq,ig,tmp,npw,npw)
	    w_0(n_pbmt+1:n_pbt,n_pbmt+1:n_pbt)=tmp
	    deallocate(tmp)
	  endif
c ------- Obtaining  W_im from W_mi -------------------------------
	  do j=1,n_pbmt
	    do i=n_pbmt+1,n_pbt
	      w_0(i,j)=conjg(w_0(j,i))
	    enddo
	  enddo
      end
