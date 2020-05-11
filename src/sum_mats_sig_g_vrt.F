      real*8 function sum_mats_sig_g_vrt(n,n0,sig,g)
c - Performs the summation of Re{Sigma(w)Gc(w)} over Matsubara frequences -
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0
	real*8, intent(in) :: sig(n0,n0,2,ndim3_omega)
	complex*16, intent(in) :: g(n0,n0,0:n_omega)
	integer :: i_omega,j,i,ind_omega
	real*8 :: de,e
	complex*16, allocatable :: tmp(:,:)
	allocate(tmp(n,n))
	e=0.d0
	do ind_omega=1,ndim3_omega
	  i_omega=me_t*ndim3_omega+ind_omega-1
        call ferm_unpack_omega(tmp,sig(1,1,1,ind_omega),n,n0,n0)
	  de=0.d0
	  do j=1,n
	    do i=1,n
	      de=de+tmp(i,j)*g(j,i,i_omega)
	    enddo
	  enddo
        e=e+et_sum_sc_gc(i_omega)*de
      enddo
	if(nproc_t/=1) call dgop(e,1,'  +',comm_t)
	sum_mats_sig_g_vrt=e
	deallocate(tmp)
      end