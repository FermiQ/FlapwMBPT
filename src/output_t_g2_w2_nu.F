      subroutine output_t_g2_w2_nu(key,t_nu,nom,i_nu)
c     key = 1 for Q
c     key = 2 for T      
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: key,nom,i_nu
	complex*16, intent(in) :: t_nu(nwp_full,nwp_full,nom,nqdiv)
	integer :: i_len,i_nu1,k
	i_len=len_trim(allfile)
	if(maswrk) then
	  if(key==1) then
	    open(3,file=allfile(1:i_len)//'_Q_bar_vertex_nu1_Re')
	    open(4,file=allfile(1:i_len)//'_Q_bar_vertex_nu1_Im')
	  else if(key==2) then
	    open(3,file=allfile(1:i_len)//'_T_vertex_nu1_Re')
	    open(4,file=allfile(1:i_len)//'_T_vertex_nu1_Im')
	  endif
	  do i_nu1=1,nom
 	    write(3,'(f11.5,5(1x,e12.5))')
     &	    w_nu_adapt_nu(i_nu1,i_nu)*evolt/2,
     &	    (real(t_nu(1,1,i_nu1,k)),k=1,min(nqdiv,5))
 	    write(4,'(f11.5,5(1x,e12.5))')
     &	    w_nu_adapt_nu(i_nu1,i_nu)*evolt/2,
     &	    (imag(t_nu(1,1,i_nu1,k)),k=1,min(nqdiv,5))
        enddo
	  close(3)
	  close(4)
	endif
      end
      
      
      subroutine output_t_g2_w1_nu(key,t_nu,nom,i_nu,np)
c     key = 1 for Q
c     key = 2 for T      
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: key,nom,i_nu,np
	complex*16, intent(in) :: t_nu(nwpmax,nwpmax,nom)
	integer :: i_len,i_nu1,k
	i_len=len_trim(allfile)
	if(maswrk) then
	  if(key==1) then
	    open(3,file=allfile(1:i_len)//'_Q_bar_vertex_nu1_Re')
	    open(4,file=allfile(1:i_len)//'_Q_bar_vertex_nu1_Im')
	  else if(key==2) then
	    open(3,file=allfile(1:i_len)//'_T_vertex_nu1_Re')
	    open(4,file=allfile(1:i_len)//'_T_vertex_nu1_Im')
	  endif
	  do i_nu1=1,nom
 	    write(3,'(f11.5,5(1x,e12.5))')
     &	    w_nu_adapt_nu(i_nu1,i_nu)*evolt/2,
     &	    (real(t_nu(k,k,i_nu1)),k=1,min(np,5))
 	    write(4,'(f11.5,5(1x,e12.5))')
     &	    w_nu_adapt_nu(i_nu1,i_nu)*evolt/2,
     &	    (imag(t_nu(k,k,i_nu1)),k=1,min(np,5))
        enddo
	  close(3)
	  close(4)
	endif
      end