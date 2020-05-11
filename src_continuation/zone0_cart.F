	  subroutine zone0_cart(q,bas,q1)
	  use models_mod
	  use solid_mod
	  implicit none
	  real*8, intent(in) :: q(3),bas(3,3)
	  real*8, intent(out) :: q1(3)
	  integer :: i,j,k
	  real*8 :: d0,d,v(3),v1(3)
	  d=100000.d0
	  v1=q
	  do k=-3,3
        do j=-3,3
	      do i=-3,3
	        v=q+i*bas(:,1)+j*bas(:,2)+k*bas(:,3)
	        d0=dot_product(v,v)
	        if(d0<d) then
	          d=d0
	          v1=v
	        endif
	      enddo
	    enddo
	  enddo
	  q1=v1
	  end
