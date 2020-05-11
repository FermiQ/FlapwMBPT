      subroutine sigx_k_solid
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
      implicit none
	  integer :: ispin,k,ind_k
	  sigx_solid_k=0.d0
      allocate(pv_x(maxel,maxel,natom))
      do ispin=1,nspin
        pv_x=(0.d0,0.d0)
        call get_pc(ispin)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          call sigma_loc_fin(ind_k,k,ispin)
          call sigma_x_div(ind_k,k,ispin)
        enddo
      enddo
      deallocate(pv_x)
      call timel('******* get_pc finished ************')
      call sigma_x_rs
      if(key_hybrid==1) sigx_solid_k=alpha_hbr*sigx_solid_k
      end
