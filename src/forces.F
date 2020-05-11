      subroutine forces
	use atom_mod
	use manager_mod
      use models_mod
	use parallel_mod
	use units_mod
	use solid_mod
      use vertex_mod
      implicit none
      real*8, allocatable :: f(:,:)
      if(ubi/='dft') return
      allocate(f(3,nsort))
      f=0.d0
      call forces_hf(f)
      call forces_core(f)
      call forces_gradro(f)
      call forces_bnd(f)
      deallocate(f)
      call timel('***** FORCES finished **************')
      end