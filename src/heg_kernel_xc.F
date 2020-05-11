      subroutine heg_kernel_xc
	use atom_mod
	use etot_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      call heg_alda(ilda,chi_cmp)
      end