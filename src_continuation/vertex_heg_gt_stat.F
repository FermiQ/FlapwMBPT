      subroutine vertex_heg_gt_stat(ispin,lambda_dyn_tau,t_tau)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      complex*16, intent(in) :: t_tau(nrr_div,ndim3_tau)
      complex*16, intent(inout) :: lambda_dyn_tau(nrr_div,2,ndim3_tau)
      integer :: i_tau,ir,ir1,ind_tau,it
      real*8 :: v(3)
      do ir=1,nrr_div
	  v=-r_mesh_heg(:,ir)
        call zone1_number(v,q0b0,nr_div,ir1)
        do ind_tau=1,ndim3_tau
	    i_tau=me_t*ndim3_tau+ind_tau-1
	    do it=1,2
	      if(it==2) i_tau=n_tau-i_tau
            lambda_dyn_tau(ir,it,ind_tau)=lambda_dyn_tau(ir,it,ind_tau)
     &                                   +g_r_tau_heg(ir,i_tau,ispin)
     &                                   *t_tau(ir,ind_tau)
          enddo
        enddo    !! over ind_tau
      enddo   !! over ir
      end