      subroutine sum_mats_sig_gc_heg(e_sig_gc,sig_omega,gc_omega,
     &                               sig_omega_as,gc_omega_as)
c -------- Performs the summation of Re{Sigma(w)Gc(w)} ---------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(in) :: sig_omega_as(2),gc_omega_as(2)
	complex*16, intent(in) :: sig_omega(ndim3_omega),
     &	                        gc_omega(ndim3_omega)
	real*8, intent(out) :: e_sig_gc
      integer :: i_omega,ind_omega
	e_sig_gc=0.d0
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        e_sig_gc=e_sig_gc
     &          +et_sum_sc_gc(i_omega)
     &          *real(sig_omega(ind_omega)*gc_omega(ind_omega))
      enddo
	if(nproc_t/=1) call dgop(e_sig_gc,1,'  +',comm_t)
c -------- Asymptotic correction is being done here -----------------
      e_sig_gc=e_sig_gc+sig_omega_as(1)*gc_omega_as(1)
     &                 *sums_omega_asymp(0,6,1)
     &                 -sig_omega_as(2)*gc_omega_as(2)
     &                 *sums_omega_asymp(0,4,1)
	e_sig_gc=e_sig_gc*2.d0
      end