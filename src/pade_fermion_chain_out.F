      subroutine pade_fermion_all_out(n,n0,ap,om,ff)
	  use manager_mod
	  use units_mod
      implicit none
      integer, intent(in) :: n,n0
      complex*16, intent(in) :: ap(0:n_omega,n0,n0),om
	  complex*16, intent(out) :: ff(n0,n0)
	  integer :: j,ii,jj
	  real*8 :: pi2b
	  complex*16 :: c1
	  complex*16, allocatable :: zi(:)
	  pi2b=(pi+pi)/betta_t
	  allocate(zi(0:n_omega))
	  zi=w_omega(0:n_omega)*(0.d0,1.d0)
c --------------------------------------------------------------
	  do jj=1,n
	    do ii=1,n
          c1=(1.d0,0.d0)
          do j=n_omega,1,-1
            c1=(1.d0,0.d0)+(om-zi(j-1))*ap(j,ii,jj)/c1
          enddo
          ff(ii,jj)=ap(0,ii,jj)/c1
        enddo
      enddo
	  deallocate(zi)
      end
