      subroutine aux_vrt_add(fif)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      real*8, intent(out) :: fif(nrel,maxel_red**2,n_pbmtm_red,nsort,
     &                           nspin)
      integer :: ispin,i,i1,isort,ij,irl
      real*8, allocatable :: fifim(:,:,:,:,:)
      fif=0.d0
c ---------------------------------------------------------------            
      allocate(fifim(nrel,maxel_red,maxel_red,n_pbmtm_red,nsort))
      do ispin=1,nspin
        call fi0_full_red(fifim,ispin,ispin)
        do isort=1,nsort
          ij=0
          do i1=1,lfunm_red(isort)
            do i=1,lfunm_red(isort)
              ij=ij+1
              do irl=1,nrel
                fif(irl,ij,:,isort,ispin)=fifim(irl,i,i1,:,isort)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(fifim)
      call timel('**** AUX_VRT_ADD finished **********')
      end
