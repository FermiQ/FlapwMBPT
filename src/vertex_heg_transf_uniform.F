      subroutine vertex_heg_transf_uniform(lambda_dyn_tau,lambda_dyn)
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: lambda_dyn_tau(nrr_div_c,2,ndim3_tau)
      complex*16, intent(out) :: lambda_dyn(nrr_div_c,0:n_omega)
      integer :: ind_tau,it
      complex*16, allocatable :: ck(:,:,:)
      allocate(ck(nrr_div_c,2,ndim3_tau))
c -------- Transform Vrt(R;tau) ---- > Vrt(K;tau) ----------
      ck=lambda_dyn_tau
      do ind_tau=1,ndim3_tau
        do it=1,2
          call from_rr_to_q_heg_cube_red(ck(1,it,ind_tau))
        enddo
      enddo  !! over i_tau
c -------- Transform Vrt(K;tau) ---- > Vrt(K;omega) ----------
      call tau_to_omega_vrt(lambda_dyn,ck,nrr_div_c)
      deallocate(ck)
      end