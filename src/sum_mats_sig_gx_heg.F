      subroutine sum_mats_sig_gx_heg(e_sig_gx,sig_omega,gx,ek,
     &                               sig_omega_as)
c -------- Performs the summation of Re{Sigma(w)Gx(w)} ---------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(in) :: ek,sig_omega_as(2)
	complex*16, intent(in) :: sig_omega(ndim3_omega),gx(ndim3_omega)
	real*8, intent(out) :: e_sig_gx
      integer :: i_omega,ind_omega
	e_sig_gx=0.d0
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        e_sig_gx=e_sig_gx+et_sum_sc_gx(i_omega)
     &                   *real(sig_omega(ind_omega)*gx(ind_omega))
      enddo
	if(nproc_t/=1) call dgop(e_sig_gx,1,'  +',comm_t)
c -------- Asymptotic correction is being done here -----------------
      e_sig_gx=e_sig_gx+sig_omega_as(1)*ek
     &                 *(sums_omega_asymp(0,4,1)
     &                  +ek*ek*sums_omega_asymp(0,6,1))
     &                 +sig_omega_as(2)
     &                 *(sums_omega_asymp(0,2,1)
     &                  -ek*ek*sums_omega_asymp(0,4,1))
	e_sig_gx=e_sig_gx*2.d0
      end