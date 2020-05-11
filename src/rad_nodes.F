      subroutine rad_nodes(phi,n,nod)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: phi(0:n)
      integer, intent(out) :: nod
	integer :: ir,i0,i1
      real*8 :: mx
      mx=maxval(abs(phi))*0.001d0
      do ir=0,n
        if(abs(phi(ir))>mx) then
          i0=ir
          goto 1
        endif
      enddo
1     do ir=n,0,-1
        if(abs(phi(ir))>mx) then
          i1=ir
          goto 2
        endif
      enddo
2	nod=0
	do ir=i0,min(i1,n-1)
	  if(phi(ir)*phi(ir+1)<0.d0) nod=nod+1
	enddo
      end