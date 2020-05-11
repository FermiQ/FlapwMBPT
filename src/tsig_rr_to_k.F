      subroutine tsig_rr_to_k(n1,n2,t_rs)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: n1,n2
      complex*16, intent(inout) :: t_rs(n1,n2,nqdiv_c)
      integer :: ir,ind,i,j,k,ij
      complex*16, allocatable :: kr(:,:)
      allocate(kr(nqdiv_c,n1*n2))
c ----- Lambda_K_TAU from Lambda_R_TAU ----------------------------
	  do ir=1,nqdiv_c
	    ij=0
        do j=1,n2
          do i=1,n1
            ij=ij+1
            kr(ir,ij)=t_rs(i,j,ir)
          enddo
        enddo
      enddo
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n1*n2,kr,-1)
	  do k=1,nqdiv_c
        call zone1_number(pnt_c(1,k),rb0_c,ndiv_c,ind)
        ij=0
        do j=1,n2
          do i=1,n1
            ij=ij+1
            t_rs(i,j,k)=kr(ind,ij)
          enddo
        enddo
      enddo
      deallocate(kr)
      end
