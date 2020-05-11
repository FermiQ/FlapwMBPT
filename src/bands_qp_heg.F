      subroutine bands_qp_heg1(x)
	use atom_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
      implicit none
      real*8, intent(in) :: x
      integer :: j,i,i_omega,ispin
	complex*16 :: a,b
	complex*16, allocatable :: tmp(:),tmp0(:)
	allocate(tmp(2))
	allocate(tmp0(0:n_omega))
c ------------------------------------------------------------
      do ispin=1,nspin
	  do i=0,n_k_heg
	    tmp=(0.d0,0.d0)
	    tmp0=(0.d0,0.d0)
	    do i_omega=0,n_omega
	      tmp0(i_omega)=sigc_k_omega_heg(i,i_omega,ispin)
	      if(i_omega<=1) then
	        j=i_omega+1
	        tmp(j)=sigc_k_omega_heg(i,i_omega,ispin)
            endif
          enddo   !! over i_omega
	    call sigma_0_z_ren(1,1,tmp,tmp0)
	    sig_0_heg(i,ispin)=tmp(1)
	    z_ren_heg(i,ispin)=tmp(2)
          b=z_ren_heg(i,ispin)
c  -------------- b is inverse Z factor now ----------------------
	    a=sig_0_heg(i,ispin)+e_star_x(i,ispin)-x
          e_star_qp(i,ispin)=a/b+x
        enddo
      enddo
	deallocate(tmp,tmp0)
      end