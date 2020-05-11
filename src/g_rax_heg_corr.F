      subroutine g_rax_heg_corr
	use atom_mod
	use etot_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: i,i_omega,ispin
	real*8 :: deltae,w
	complex*16 :: c1,sig,g_om
      complex*16, allocatable :: achain(:)
      allocate(achain(0:n_omega))
      deltae=(emaxdos-emindos)/ndos
      do ispin=1,nspin
	  do i=0,n_k_heg
          call pade_fermion_all(sigc_k_omega_heg(i,:,ispin),1,1,achain)
	    do i_omega=0,ndos
            w=emindos+i_omega*deltae
	      c1=dcmplx(w,0.d0)
		  call pade_fermion_all_out(1,1,achain,c1,sig)
		  sig_real_heg(i,i_omega,ispin)=sig
		  g_om=(0.d0,0.d0)
	      c1=dcmplx(w,e_small)
		  g_om=(1.d0,0.d0)/(c1+chem_pot-e_star_x(i,ispin)-sig)
            if(ubi/=' qp') g_rax_heg(i_omega,i,ispin)=g_om
          enddo
        enddo
      enddo
	deallocate(achain)
      end