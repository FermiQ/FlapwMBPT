      subroutine k_k_to_rr_all(kf,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: n
      complex*16, intent(inout) :: kf(n,nqdiv_c)
      integer :: ir,ind,i,k
      complex*16, allocatable :: kr(:,:)
      allocate(kr(nqdiv_c,n))
      do k=1,nqdiv_c
        call zone1_number(pnt_c(1,k),rb0_c,ndiv_c,ind)
        do i=1,n
            kr(ind,i)=kf(i,k)
        enddo
      enddo
      call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n,kr,1)
      do ir=1,nqdiv_c
        do i=1,n
          kf(i,ir)=kr(ir,i)/nqdiv_c
        enddo
      enddo
      deallocate(kr)
      end
