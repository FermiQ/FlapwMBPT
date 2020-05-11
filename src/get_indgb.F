      subroutine get_indgb
	use manager_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer :: nbndt,k
	real*8 :: cut
	do k=1,npnt
	  call g_inside_cutoff(nplwbas,pnt(1,k),cutoff,.true.,.true.,
     &	                   indgb(1,k),nbndt,nbndtm,nrdiv)
	enddo
	do k=1,num_k_plot
	  call g_inside_cutoff(nplwbas,k_plot(1,k),cutoff,.true.,.true.,
     &	                   indgb_plot(1,k),nbndt,nbndtm,nrdiv)
	enddo
      if(mesh_k_c) then
	  allocate(indgb_red(nbasmpw_red/nrel,nqdiv))
	  cut=min(cutoff,cut_bnd_red)
	  do k=1,npnt
	    call g_inside_cutoff(nplwbas,pnt(1,k),cut,.true.,.true.,
     &	                     indgb_red(1,k),nbndt,nbndtm,nrdiv_red)
	  enddo
      endif
      end