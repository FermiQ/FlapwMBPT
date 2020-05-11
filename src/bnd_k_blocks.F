      subroutine bnd_k_blocks
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use units_mod
	  use solid_mod
      implicit none
	  integer :: i,i1,i0,ndeg,ispin,k,ibl,irr
	  real*8 :: eps
	  eps=1.d-10
	  if(allocated(n_block_k)) deallocate(n_block_k)
	  allocate(n_block_k(npnt,nspin))
	  k_irrep_max=0
	  do ispin=1,nspin
	    do k=1,npnt
	      ibl=0
	      irr=0
	      i1=0
1	      i0=i1+1
	      ndeg=1
	      do i=i0+1,n_bnd(k,ispin)
            if(abs(e_bnd(i,k,ispin)-e_bnd(i0,k,ispin))<eps) ndeg=ndeg+1
	      enddo
	      i1=i0+ndeg-1
	      ibl=ibl+1
	      irr=irr+ndeg*ndeg
          if(i1<n_bnd(k,ispin)) goto 1
          n_block_k(k,ispin)=ibl
          k_irrep_max=max(k_irrep_max,irr)
        enddo
      enddo
      n_block_k_max=maxval(n_block_k)
	  if(allocated(ind_block_k)) deallocate(ind_block_k)
	  if(allocated(last_block_k)) deallocate(last_block_k)
	  if(allocated(entry_block_k)) deallocate(entry_block_k)
      allocate(ind_block_k(nbndf,npnt,nspin))
      allocate(last_block_k(n_block_k_max,npnt,nspin))
      allocate(entry_block_k(n_block_k_max,npnt,nspin))
      ndeg_max=0
	  do ispin=1,nspin
	    do k=1,npnt
	      ibl=0
	      irr=0
	      i1=0
2	      i0=i1+1
	      ndeg=1
	      do i=i0+1,n_bnd(k,ispin)
	        if(abs(e_bnd(i,k,ispin)-e_bnd(i0,k,ispin))<eps)
     &	      ndeg=ndeg+1
	      enddo
	      ndeg_max=max(ndeg_max,ndeg)
	      i1=i0+ndeg-1
	      ibl=ibl+1
	      ind_block_k(i0:i1,k,ispin)=ibl
	      last_block_k(ibl,k,ispin)=i1
	      entry_block_k(ibl,k,ispin)=irr+1
	      irr=irr+ndeg*ndeg
          if(i1<n_bnd(k,ispin)) goto 2
        enddo
      enddo
      end
