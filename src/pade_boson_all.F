      subroutine pade_boson_all(ff,wrx,nrax,xir,key)
c     key = 1 - for functions with 1/nu asymptotic
c     key = 2 - for functions with 1/nu^2 asymptotic      
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: nrax,key
      real*8, intent(in) :: wrx(0:nrax)
	  complex*16, intent(in) :: ff(0:n_nu)
	  complex*16, intent(out) :: xir(0:nrax)
	  integer :: n,i,j
      real*8 :: pi2b
	  complex*16 :: cnu,c1
	  complex*16, allocatable :: a(:),f(:),g(:),zi(:)
      pi2b=(pi+pi)/betta_t
	  n=n_nu
	  allocate(a(0:n))
	  allocate(f(0:n))
	  allocate(g(0:n))
	  allocate(zi(0:n))
      do i=0,n
        zi(i)=dcmplx(0.d0,w_nu(i))**key
      enddo
	  a(0)=ff(0)
	  f=ff(0:n)
	  do j=1,n
	    do i=j,n
	      g(i)=(a(j-1)-f(i))/(zi(i)-zi(j-1))/f(i)
	    enddo
	    a(j)=g(j)
	    f=g
	  enddo
c --------------------------------------------------------------
      do i=0,nrax
        cnu=dcmplx(wrx(i),e_small_bos)**key
        c1=(1.d0,0.d0)
        do j=n,1,-1
          c1=(1.d0,0.d0)+(cnu-zi(j-1))*a(j)/c1
        enddo
        xir(i)=a(0)/c1
      enddo
	  deallocate(a,f,g,zi)
      end
