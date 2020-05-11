      subroutine output_p_static_q(p00,pzz,pw0,pwz)
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      real*8, intent(in) :: pw0(nc_line),pwz(nc_line)
	  complex*16, intent(in) :: p00(n_pbtot_red,n_pbtot_red,nc_line),
     &                          pzz(n_pbtot_red,n_pbtot_red,nc_line)
	  integer :: i_len,k,ns,j,i,k0,ka
	  real*8 :: q,pi2a
	  real*8, allocatable :: p(:,:)
	  pi2a=(pi+pi)/par
	  allocate(p(nc_line,2))
	  p=0.d0
	  do k=1,nc_line
	    k0=kline_in_npnt_c(k)
        ka=k_a_from_c(k0)
	    ns=n_pbmt_red+nplwgw_red(ka)
	    do j=1,ns
	      do i=1,ns
	        p(k,1)=p(k,1)+conjg(pw_pb_red(i,ka))*p00(i,j,k)
     &                   *pw_pb_red(j,ka)
	        p(k,2)=p(k,2)+conjg(pw_pb_red(i,ka))*pzz(i,j,k)
     &                   *pw_pb_red(j,ka)
          enddo
        enddo
      enddo
	  i_len=len_trim(allfile)
	  if(maswrk) then
	    if(ubi==' hf') open(3,file=allfile(1:i_len)//
     &	                    '_P_Static_vertex_q.hf')
	    if(ubi==' gw') open(3,file=allfile(1:i_len)//
     &	                    '_P_Static_vertex_q.gw')
	    if(ubi==' qp') open(3,file=allfile(1:i_len)//
     &	                    '_P_Static_vertex_q.qp')
	    if(ubi=='psi') open(3,file=allfile(1:i_len)//
     &	                    '_P_Static_vertex_q.psi')
	    if(ubi=='bsp') open(3,file=allfile(1:i_len)//
     &	                    '_P_Static_vertex_q.bsp')
	    write(3,*)' d_P_00 and d_P_ZZ from static Chi-calculations'
        write(3,*)' dP00_PB  dP00_PW  dPZZ_PB  dPZZ_PW in (1/Ry/a.u.v.)'
	    do k=1,nc_line
	      k0=kline_in_npnt_c(k)
	      q=sqrt(dot_product(pnt_c(:,k0),pnt_c(:,k0)))*pi2a
 	      write(3,'(f11.5,4(1x,e12.5))')q,p(k,1),pw0(k),p(k,2),pwz(k)
	    enddo
	    close(3)
	  endif
	  deallocate(p)
      end
