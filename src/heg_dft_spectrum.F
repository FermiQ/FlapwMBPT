      subroutine heg_dft_spectrum
	  use heg_mod
	  use etot_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer :: i,ispin
      real*8 :: vxcup
      call EXCVXC(ro_spin(2),ro_spin(1),ilda,EXC_heg(2),EXC_heg(1),
     &            vxc_heg(2),vxc_heg(1))
      CALL LSDEXCH(ro_el_gas,ex_heg,vxcup)
c ------------------------------------------------------------------
      do ispin=1,nspin
        do i=0,n_k_heg
          e_star_x(i,ispin)=e_star_free(i,ispin)+vxc_heg(ispin)
        enddo
      enddo
      end
