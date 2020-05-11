      subroutine pade_fermion_all(ff,n,n0,ap)
	  use manager_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0
	  complex*16, intent(in) :: ff(n0,n0,0:n_omega)
	  complex*16, intent(out) :: ap(0:n_omega,n0,n0)
	  integer :: i,j,ii,jj
	  real*8 :: sm,pib
	  complex*16, allocatable :: f(:),g(:),zi(:)
      pib=pi/betta_t
	  allocate(f(0:n_omega))
	  allocate(g(0:n_omega))
	  allocate(zi(0:n_omega))
      do i=0,n_omega
	    zi(i)=w_omega(i)*(0.d0,1.d0)
      enddo
	  ap=(0.d0,0.d0)
	  do jj=1,n
	    do ii=1,n
	      sm=sum(abs(ff(ii,jj,:)))
	      if(sm<1.d-12) cycle
	      ap(0,ii,jj)=ff(ii,jj,0)
	      f=ff(ii,jj,0:n_omega)
	      do j=1,n_omega
	        do i=j,n_omega
	          g(i)=(ap(j-1,ii,jj)-f(i))/(zi(i)-zi(j-1))/f(i)
	        enddo
	        ap(j,ii,jj)=g(j)
	        f=g
	      enddo
	    enddo
	  enddo
	  deallocate(f,g,zi)
      end
