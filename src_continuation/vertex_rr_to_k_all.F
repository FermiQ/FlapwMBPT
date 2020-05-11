      subroutine vertex_rr_to_k_all(vx,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: n
      complex*16, intent(inout) :: vx(n,nqdiv_c)
      integer :: ir,ind,i,k
      complex*16, allocatable :: kr(:,:)
      allocate(kr(nqdiv_c,n))
	  do ir=1,nqdiv_c
	    do i=1,n
          kr(ir,i)=vx(i,ir)
        enddo
      enddo
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n,kr,-1)
	  do k=1,nqdiv_c
        call zone1_number(pnt_c(1,k),rb0_c,ndiv_c,ind)
        do i=1,n
          vx(i,k)=kr(ind,i)
        enddo
      enddo
      deallocate(kr)
      end
