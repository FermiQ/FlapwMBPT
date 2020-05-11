      subroutine output_p_vertex_red(title,p_line)
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      character*8, intent(in) :: title
	  real*8, intent(in) :: p_line(n_pbtot_red,n_pbtot_red,ndim3_nu,
     & 	                           n_line)
	  integer :: i_len,np,k,i_nu,j,i,ind_nu,k0,nl
	  real*8, allocatable :: p(:,:)
	  complex*16, allocatable :: pc(:,:)
	  nl=min(5,n_line)
	  allocate(p(0:n_nu,nl))
	  allocate(pc(n_pbtot_red,n_pbtot_red))
	  p=0.d0
      do ind_nu=1,ndim3_nu
	    i_nu=me_t*ndim3_nu+ind_nu-1
	    do k=1,nl
	      k0=kline_in_npnt(k)
	      np=n_pbmt_red+nplwgw_red(k0)
	      call unpack_hermit(pc,p_line(1,1,ind_nu,k),np,n_pbtot_red,
     &	                     n_pbtot_red,(0.d0,0.d0),(1.d0,0.d0))
	      do j=1,np
	        do i=1,np
	          p(i_nu,k)=p(i_nu,k)+conjg(pw_pb_red(i,k0))*pc(i,j)
     &                                 *pw_pb_red(j,k0)
            enddo
          enddo
        enddo
      enddo
	  if(nproc_t/=1) call dgop(p,(n_nu+1)*nl,'  +',comm_t)
	  p_q0_vrt=p(0:1,1)
	  i_len=len_trim(allfile)
	  if(maswrk) then
	    open(3,file=allfile(1:i_len)//'_'//title)
	    do i_nu=0,n_nu
 	      write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                                (p(i_nu,k),k=1,nl)
	    enddo
	    close(3)
	  endif
	  deallocate(p,pc)
      end
