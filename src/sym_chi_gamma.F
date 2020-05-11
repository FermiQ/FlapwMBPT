      subroutine sym_chi_gamma(ig,w_0,n)
	use atom_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ig,n
	complex*16, intent(inout) :: w_0(n,n)
	integer :: j,i
	complex*16, allocatable :: tmp(:,:)
	call sym_w_mt(ig,pnt,w_0,n)
	allocate(tmp(n_pbmt,nplw_gw))
	tmp(:,1:nplwgw(1))=w_0(1:n_pbmt,n_pbmt+1:n)
	call sym_w_mi(1,ig,tmp)
	w_0(1:n_pbmt,n_pbmt+1:n)=tmp(:,1:nplwgw(1))
	deallocate(tmp)
	allocate(tmp(nplwgw(1),nplwgw(1)))
	tmp=w_0(n_pbmt+1:n,n_pbmt+1:n)
	call sym_w_ii(1,ig,tmp,nplwgw(1),nplwgw(1))
	w_0(n_pbmt+1:n,n_pbmt+1:n)=tmp
	deallocate(tmp)
c ------- Obtaining  W_im from W_mi -------------------------------
	do j=1,n_pbmt
	  do i=n_pbmt+1,n
	    w_0(i,j)=conjg(w_0(j,i))
	  enddo
	enddo
      end
      
      
      
      subroutine sym_chi_q(ig,w_0,n,iq)
	  use atom_mod
	  use solid_mod
      implicit none
      integer, intent(in) :: ig,n,iq
	  complex*16, intent(inout) :: w_0(n,n)
	  integer :: j,i,npw
	  complex*16, allocatable :: tmp(:,:)
      npw=n-n_pbmt
	  call sym_w_mt(ig,pnt(1,iq),w_0,n)
	  allocate(tmp(n_pbmt,nplw_gw))
	  tmp(:,1:npw)=w_0(1:n_pbmt,n_pbmt+1:n)
	  call sym_w_mi(iq,ig,tmp)
	  w_0(1:n_pbmt,n_pbmt+1:n)=tmp(:,1:npw)
	  deallocate(tmp)
	  allocate(tmp(npw,npw))
	  tmp=w_0(n_pbmt+1:n,n_pbmt+1:n)
	  call sym_w_ii(iq,ig,tmp,npw,npw)
	  w_0(n_pbmt+1:n,n_pbmt+1:n)=tmp
	  deallocate(tmp)
c ------- Obtaining  W_im from W_mi -------------------------------
	  do j=1,n_pbmt
	    do i=n_pbmt+1,n
	      w_0(i,j)=conjg(w_0(j,i))
	    enddo
	  enddo
      end
