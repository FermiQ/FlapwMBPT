      subroutine sym_w_red(iq,iq0,w_0)
	  use atom_mod
	  use solid_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,iq0
	  complex*16, intent(inout) :: w_0(n_pbtot_red,n_pbtot_red)
	  integer :: j,i,ig,ka0
	  complex*16, allocatable :: tmp(:,:)
c -------- Copy+Symmetrization of W into another k-point ---------------
	  ig=k_group_c(iq)
      if(iq/=iq0) then
	    call sym_w_mt_red(ig,pnt_c(1,iq),w_0,n_pbtot_red)
	    allocate(tmp(n_pbmt_red,nplw_gw_red))
	    tmp=w_0(1:n_pbmt_red,n_pbmt_red+1:n_pbtot_red)
	    call sym_w_mi_red(iq,ig,tmp)
	    w_0(1:n_pbmt_red,n_pbmt_red+1:n_pbtot_red)=tmp
	    deallocate(tmp)
	    allocate(tmp(nplw_gw_red,nplw_gw_red))
	    tmp=w_0(n_pbmt_red+1:n_pbtot_red,n_pbmt_red+1:n_pbtot_red)
        ka0=k_a_from_c(iq0)
	    call sym_w_ii_red(iq,ig,tmp,nplwgw_red(ka0),nplw_gw_red)
	    w_0(n_pbmt_red+1:n_pbtot_red,n_pbmt_red+1:n_pbtot_red)=tmp
	    deallocate(tmp)
	  endif
c ------- Obtaining  W_im from W_mi -------------------------------
	  do j=1,n_pbmt_red
	    do i=n_pbmt_red+1,n_pbtot_red
	      w_0(i,j)=conjg(w_0(j,i))
	    enddo
	  enddo
      end
