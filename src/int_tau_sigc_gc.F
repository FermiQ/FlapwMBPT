      subroutine int_tau_sigc_gc(e,f1_tau,f2_tau,n,n0,n1)
c ------ Performs the integral of 1/2 F1(tau)*F2(beta-tau) ---------
c --------------------  TAU-functions are given in A;B forms ----
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f1_tau(n0,n0,2,ndim3_tau),
     &	                    f2_tau(n1,n1,2,ndim3_tau)
	  real*8, intent(out) :: e
	  integer :: ind_tau,i0,i
	  real*8 :: conv_hh
	  real*8, allocatable :: f(:)
	  allocate(f(ndim3_tau))
c ------------- AA and BB parts --------------------------------------
      do ind_tau=1,ndim3_tau
        f(ind_tau)=conv_hh(f1_tau(1,1,1,ind_tau),f2_tau(1,1,1,ind_tau),
     &                     n,n0,n1)
     &            -conv_hh(f1_tau(1,1,2,ind_tau),f2_tau(1,1,2,ind_tau),
     &                     n,n0,n1)
      enddo
      i0=me_t*ndim3_tau
      i=i0+ndim3_tau-1
      e=dot_product(f,ab_nu_from_tau(i0:i,0,1))/4.d0
	  deallocate(f)
      end
