      subroutine bands_vertex
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ispin,k,ibnd,ib,mi,ma
      integer, allocatable :: ind_bands(:,:,:)
      allocate(ind_bands(nbndf,npnt,nspin))
      ind_bands=0
      if(allocated(n_low_bnd)) deallocate(n_low_bnd)
      allocate(n_low_bnd(npnt,nspin_0))
      n_low_bnd=0
      if(maswrk) write(iun,*)' Bands included in VERTEX calculation'
c ----------------------------------------------------------
      do ispin=1,nspin
        do k=1,npnt
          call bands_near_ef(nbndf_bnd_0(1),nbndf_bnd_0(2),
     &	                   n_bnd(k,ispin),e_bnd(1,k,ispin),mi,ma)
          n_low_bnd(k,ispin)=ma-mi+1
          do ib=mi,ma
            ibnd=ib-mi+1
            ind_bands(ibnd,k,ispin)=ib
          enddo
          ib=n_low_bnd(k,ispin)
          if(maswrk) write(iun,'(a10,2i4,a19)')
     &		        ' Ispin k ',ispin,k,' Low energy Bands :'
          if(maswrk) write(iun,'(17i4)')(ind_bands(ibnd,k,ispin),
     &	                               ibnd=1,ib)
        enddo
      enddo
      nbndf_bnd=maxval(n_low_bnd)
      if(allocated(ind_bands_bnd)) deallocate(ind_bands_bnd)
      allocate(ind_bands_bnd(nbndf_bnd,npnt,nspin_0))
      do ispin=1,nspin
        do k=1,npnt
          do ibnd=1,n_low_bnd(k,ispin)
            ind_bands_bnd(ibnd,k,ispin)=ind_bands(ibnd,k,ispin)
          enddo
	    call size_shift_par(n_low_bnd(k,ispin),nproc_b,
     &                        ndim4_bnd(1,k,ispin),
     &                        n4_mpi_bnd(1,k,ispin))
        enddo
      enddo
      deallocate(ind_bands)
      end