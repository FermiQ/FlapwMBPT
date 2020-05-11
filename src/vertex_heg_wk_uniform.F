      subroutine vertex_heg_wk_uniform(k_pw,lambda_dyn_tau,lambda_stat)
	use atom_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: k_pw(nrr_div_c,2,ndim3_tau)
      real*8, intent(out) :: lambda_dyn_tau(nrr_div_c,2,ndim3_tau),
     &                       lambda_stat(nrr_div_c)
      integer :: i,i0_tau,ind_tau,it
      complex*16, allocatable :: tm(:)
      lambda_stat=0.d0
      allocate(tm(nrr_div_c))
c ----- Lambda = W*K ----------------------------
	do ind_tau=1,ndim3_tau
	  i0_tau=me_t*ndim3_tau+ind_tau-1
	  do it=1,2
          do i=1,nrr_div_c
            lambda_dyn_tau(i,it,ind_tau)=w_r_tau_heg_c(i,i0_tau)
     &                                  *k_pw(i,it,ind_tau)
          enddo
        enddo
        if(i0_tau==0) then
          do i=1,nrr_div_c
            lambda_stat(i)=v_r_heg_c(i)*k_pw(i,1,1)
          enddo
          tm=lambda_stat
          call from_rr_to_q_heg_cube_red(tm)
          lambda_stat=tm
        endif
      enddo  !! over i_tau
      deallocate(tm)
	if(nproc_t/=1) call dgop(lambda_stat,nrr_div_c,'  +',comm_t)
      end