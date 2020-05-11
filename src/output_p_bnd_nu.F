      subroutine output_p_bnd_nu(key,pol,n0)
c     key = 1 for P_0
c     key = 2 for delta_P     
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: key,n0
	  complex*16, intent(in) :: pol(n0,n0,ndim3_nu,
     &                              ndim_k_red(me_k+1))
	  integer :: i_len,np,k,ns,i_nu,j,i,ind_nu,k0,ind_k,ka
	  real*8 :: q
	  real*8, allocatable :: p(:,:)
	  np=min(5,nc_line)
	  allocate(p(0:n_nu,nc_line))
	  p=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        do ind_k=1,ndim_k_red(me_k+1)
          k0=n_mpi_k_red(me_k+1)+ind_k
	      k=k_npnt_in_line_c(k0)
	      if(k==0) cycle
          ka=k_a_from_c(k0)
	      ns=n_pbmt_red
	      if(n0>n_pbmt_red) ns=n_pbmt_red+nplwgw_red(ka)
	      do j=1,ns
	        do i=1,ns
	          p(i_nu,k)=p(i_nu,k)+conjg(pw_pb_red(i,ka))
     &	                         *pol(i,j,ind_nu,ind_k)
     &                           *pw_pb_red(j,ka)
            enddo
          enddo
        enddo
      enddo
	  if(nproc_t/=1) call dgop(p,(n_nu+1)*nc_line,'  +',comm_t)
	  if(nproc_k/=1) call dgop(p,(n_nu+1)*nc_line,'  +',comm_k)
	  i_len=len_trim(allfile)
	  if(maswrk) then
	    if(key==1) open(3,file=allfile(1:i_len)//'_P_0_vertex_nu')
	    if(key==2) open(3,file=allfile(1:i_len)//'_P_1_vertex_nu')
	    do i_nu=0,n_nu
 	      write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu),(p(i_nu,k),k=1,np)
	    enddo
	    close(3)
	    if(key==1) open(3,file=allfile(1:i_len)//'_P_0_vertex_q')
	    if(key==2) open(3,file=allfile(1:i_len)//'_P_1_vertex_q')
	    do k=1,nc_line
	      k0=kline_in_npnt_c(k)
	      q=sqrt(dot_product(pnt_c(:,k0),pnt_c(:,k0)))
 	      write(3,'(f11.5,5(1x,e12.5))')q,(p(i_nu,k),i_nu=0,4)
	    enddo
	    close(3)
	  endif
	  deallocate(p)
      end
      
      subroutine output_p_bnd_nu_0(key,pol,n0)
c     key = 1 for P_0
c     key = 2 for delta_P     
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: key,n0
	  complex*16, intent(in) :: pol(n0,n0,0:n_nu,npnt_c)
	  integer :: i_len,np,k,ns,i_nu,j,i,ind_nu,k0,ka
	  real*8 :: q
	  real*8, allocatable :: p(:,:)
	  np=min(5,nc_line)
	  allocate(p(0:n_nu,nc_line))
	  p=0.d0
      do ind_nu=1,ndim3_nu
	    i_nu=me_t*ndim3_nu+ind_nu-1
	    do k=1,nc_line
	      k0=kline_in_npnt_c(k)
          ka=k_a_from_c(k0)
	      ns=n_pbmt_red+nplwgw_red(ka)
	      do j=1,ns
	        do i=1,ns
	          p(i_nu,k)=p(i_nu,k)+conjg(pw_pb_red(i,ka))
     &	                         *pol(i,j,i_nu,k0)
     &                           *pw_pb_red(j,ka)
            enddo
          enddo
        enddo
      enddo
	  if(nproc_t/=1) call dgop(p,(n_nu+1)*nc_line,'  +',comm_t)
	  i_len=len_trim(allfile)
	  if(maswrk) then
	    if(key==1) open(3,file=allfile(1:i_len)//'_P_0_vertex_nu')
	    if(key==2) open(3,file=allfile(1:i_len)//'_P_1_vertex_nu')
	    do i_nu=0,n_nu
 	      write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                                (p(i_nu,k),k=1,np)
	    enddo
	    close(3)
	    if(key==1) open(3,file=allfile(1:i_len)//'_P_0_vertex_q')
	    if(key==2) open(3,file=allfile(1:i_len)//'_P_1_vertex_q')
	    do k=1,nc_line
	      k0=kline_in_npnt_c(k)
	      q=sqrt(dot_product(pnt_c(:,k0),pnt_c(:,k0)))
 	      write(3,'(f11.5,5(1x,e12.5))')q,(p(i_nu,k),i_nu=0,4)
	    enddo
	    close(3)
	  endif
	  deallocate(p)
      end
