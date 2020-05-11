      subroutine vertex_heg_gt(ispin,lambda_dyn_tau,t_tau)
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
      complex*16, intent(in) :: t_tau(nrr_div_c,0:n_tau/2)
      complex*16, intent(inout) :: lambda_dyn_tau(nrr_div_c,0:n_tau)
      integer :: i_tau,ir,ir1,i1_tau
      real*8 :: v(3)
      do ir=1,nrr_div_c
	  v=-r_mesh_heg_c(:,ir)
        call zone1_number(v,q0b0_c,nr_div_c,ir1)
        do i_tau=0,n_tau
          i1_tau=n_tau-i_tau
          if(i_tau<=n_tau/2) then
            lambda_dyn_tau(ir,i_tau)=lambda_dyn_tau(ir,i_tau)
     &                              +g_r_tau_heg_c(ir,i_tau,ispin)
     &                              *t_tau(ir,i_tau)
          else
            lambda_dyn_tau(ir,i_tau)=lambda_dyn_tau(ir,i_tau)
     &                              +g_r_tau_heg_c(ir,i_tau,ispin)
     &                              *conjg(t_tau(ir1,i1_tau))
          endif
        enddo    !! over i_tau
      enddo   !! over ir
      end