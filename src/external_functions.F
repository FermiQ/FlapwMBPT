      real*8 function ferm_dir_0(k)
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: k
	real*8 :: fermi_dirac
	ferm_dir_0=k*fermi_dirac(k*k-chem_pot_heg_0)
      end