      subroutine v4_for_bnd_mm(v4_r,v2_q,npb)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: npb
      complex*16, intent(in) :: v2_q(npb,npb,nqdiv_c)
	  real*8, intent(out) :: v4_r(n_pbmt_red,n_pbmt_red,nqdiv_c)
	  integer :: ir,j,i,iq,k11
	  complex*16, allocatable :: v2r(:,:,:),tmp(:,:)
c ----------------------------------------------------------------
      v4_r=0.d0
      allocate(v2r(nqdiv_c,n_pbmt_red,n_pbmt_red))
      allocate(tmp(n_pbmt_red,n_pbmt_red))	 
      v2r=(0.d0,0.d0)
      do iq=1,nqdiv_c
	    call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
	    v2r(k11,:,:)=v2_q(1:n_pbmt_red,1:n_pbmt_red,iq)
      enddo   !! over iq
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red**2,v2r,1)
	  v2r=v2r/nqdiv_c
      do ir=1,nqdiv_c
	    do j=1,n_pbmt_red
	      do i=1,n_pbmt_red
	        v4_r(i,j,ir)=v2r(ir,i,j)
	      enddo
	    enddo
      enddo   !! over ir
      deallocate(v2r,tmp)
      end
