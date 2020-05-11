      subroutine output_e_tot_vrt(e_dyn_o,e_dyn_t,e_semi_o,e_semi_t,
     &                            e_stat_o,e_stat_t,e0_dyn_o,e0_dyn_t,
     &                            e0_semi_o,e0_semi_t,e0_stat_o,
     &                            e0_stat_t)
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      real*8, intent(in) :: e_dyn_o,e_dyn_t,e_semi_o,e_semi_t,e_stat_o,
     &                      e_stat_t,e0_dyn_o,e0_dyn_t,e0_semi_o,
     &                      e0_semi_t,e0_stat_o,e0_stat_t
	  real*8 :: e0_o,e0_t,e_o,e_t
	  e0_o=e0_dyn_o+e0_semi_o+e0_stat_o
	  e0_t=e0_dyn_t+e0_semi_t+e0_stat_t
	  e_o=e_dyn_o+e_semi_o+e_stat_o
	  e_t=e_dyn_t+e_semi_t+e_stat_t
	  if(.not.maswrk) return
      write(iun,*)' '
	  write(iun,'(10x,a43)')
     &	'Internal Energy components from Vertex Part'
      write(iun,*)' '
	  write(iun,'(a10,13x,a5,21x,a3)')'Components','Omega','Tau'
	  write(iun,'(15x,a8,7x,a4,7x,a8,7x,a4)')'One-site','Full',
     &	                                     'One-site','Full'
      write(iun,*)' '
      write(iun,'(1x,a6,6x,2f12.8,2x,2f12.8)')'E_Dyn ',e0_dyn_o,
     &                                        e_dyn_o,e0_dyn_t,e_dyn_t
      write(iun,'(1x,a6,6x,2f12.8,2x,2f12.8)')'E_Semi',e0_semi_o,
     &                                     e_semi_o,e0_semi_t,e_semi_t
      write(iun,'(1x,a6,6x,2f12.8,2x,2f12.8)')'E_Stat',e0_stat_o,
     &                                     e_stat_o,e0_stat_t,e_stat_t
      write(iun,'(1x,a6,6x,2f12.8,2x,2f12.8)')'E_Full',e0_o,e_o,e0_t,
     &                                        e_t
      end
