      subroutine green_0(correlation)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
      logical, intent(in) :: correlation
	  integer :: k,ispin,len,k0,i
	  real*8 :: de,green0
c ------------- Forming G_FULL(0-) --> G_FULL_0 ----------------------
      if(.not.correlation) g_full_00=0.d0
	  if(me_t==0) then
	    do ispin=1,nspin
	      do k=1,ndim3_k(me_k+1)
	        k0=n3_mpi_k(me_k+1)+k
	        if(.not.correlation) then
	          g_full_00(:,k,ispin)=0.d0
              do i=1,n_bnd(k0,ispin)
                de=e_bnd(i,k0,ispin)-chem_pot
	            g_full_00(i,k,ispin)=-green0(de,tau_mesh(n_tau))
              enddo
            endif
	      enddo   !! over k
	    enddo  !! over ispin
	  endif
	  if(nproc_t/=1) then
	    if(.not.correlation) then
	      len=nbndf*ndim3_k(me_k+1)*nspin_0
	      if(len/=0) call dgop(g_full_00,len,'  +',comm_t)
        endif
	  endif
      end
