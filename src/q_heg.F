      subroutine q_heg(q_tau,k_tau)
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use vertex_mod
      implicit none
      complex*16, intent(in) :: k_tau(nrr_div_c,0:n_tau,nspin)
      complex*16, intent(out) :: q_tau(nrr_div_c,0:n_tau/2)
      integer :: i_tau,i1_tau,ir,ir1,ispin
      real*8 :: v(3)
c ----- Q(r;tau)[q;nu_d] ----------------------------
      q_tau=(0.d0,0.d0)
      do ir=1,nrr_div_c
	  v=-r_mesh_heg_c(:,ir)
        call zone1_number(v,q0b0_c,nr_div_c,ir1)
        do i_tau=0,n_tau/2
          i1_tau=n_tau-i_tau
          do ispin=1,nspin
            q_tau(ir,i_tau)=q_tau(ir,i_tau)
     &        -g_r_tau_heg_c(ir,i1_tau,ispin)*k_tau(ir,i_tau,ispin)
     &        -g_r_tau_heg_c(ir,i_tau,ispin)
     &        *conjg(k_tau(ir1,i1_tau,ispin))
          enddo
        enddo    !! over i_tau
      enddo   !! over ir
      if(nspin==1) q_tau=2.d0*q_tau
c ----- Q(q';tau)[q;nu] ----------------------------
      do i_tau=0,n_tau/2
        call from_rr_to_q_heg_cube_red(q_tau(1,i_tau))
      enddo  !! over i_tau
      end