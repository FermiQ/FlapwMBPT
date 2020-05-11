      subroutine output_w
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      complex*16, allocatable :: f(:,:)
	  if(.not.maswrk) return
	  if(mbpt_gw/=0) return
      allocate(f(0:n_nu,n_line))
c ------- G=G'=0 components for q-points along K-line -----------
      f=p_output
      call output_boson_q_nu(f,'_____P',n_line,'a')
      f=e_output
      if(mode_bts==0) call output_boson_q_nu(f,'_____E',n_line,'a')
      if(mode_bts==1) call output_boson_q_nu(f,'____E1',n_line,'a')
      f=w_output
      call output_boson_q_nu(f,'_____W',n_line,'a')
      deallocate(f)
      end
