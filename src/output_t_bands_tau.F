      subroutine output_t_bands_tau(key,t_pb)
c     key = 1 for Q_bar
c     key = 2 for T      
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: key
	complex*16, intent(in) :: t_pb(npb_g2_w2_max,npb_g2_w2_max,
     &	                             0:n_tau,nqdiv)
	integer :: i_len,i_tau,k
	i_len=len_trim(allfile)
	if(maswrk) then
	  if(key==1) then
	    open(3,file=allfile(1:i_len)//'_Q_bar_vertex_Re')
	    open(4,file=allfile(1:i_len)//'_Q_bar_vertex_Im')
	  else if(key==2) then
	    open(3,file=allfile(1:i_len)//'_T_vertex_Re')
	    open(4,file=allfile(1:i_len)//'_T_vertex_Im')
	  endif
	  do i_tau=0,n_tau
 	    write(3,'(f11.5,5(1x,e12.5))')tau_mesh(i_tau),
     &	      (real(t_pb(1,1,i_tau,k)),k=1,min(nqdiv,5))
 	    write(4,'(f11.5,5(1x,e12.5))')tau_mesh(i_tau),
     &	      (imag(t_pb(1,1,i_tau,k)),k=1,min(nqdiv,5))
	  enddo
	  close(3)
	  close(4)
	endif
      end