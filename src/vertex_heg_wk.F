      subroutine vertex_heg_wk(q,k_pw,lambda_dyn_tau,lambda_stat)
	use atom_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: q(*)
      complex*16, intent(in) :: k_pw(nrr_div_c,0:n_tau)
      complex*16, intent(out) :: lambda_dyn_tau(nrr_div_c,0:n_tau),
     &                           lambda_stat(nrr_div_c)
      integer :: i,i_tau,i1_tau,i0_tau,k,kq,k0,kq0
      real*8 :: v(3)
      complex*16, allocatable :: b(:)
      allocate(b(nrr_div_c))
c ----- Lambda = W*K ----------------------------
      do i_tau=0,n_tau
        i1_tau=n_tau-i_tau
        i0_tau=min(i_tau,i1_tau)
        b=k_pw(:,i_tau)
        do i=1,nrr_div_c
          lambda_dyn_tau(i,i_tau)=w_r_tau_heg_c(i,i0_tau)*b(i)
        enddo
        if(i_tau==0) then
          do i=1,nrr_div_c
            lambda_stat(i)=v_r_heg_c(i)*b(i)
          enddo
          call from_rr_to_q_heg_cube_red(lambda_stat)
        endif
      enddo  !! over i_tau
c -------- Lambda_stat = A(k)+B(k) = B*(-k+q) + B(k) -------------
      b=lambda_stat
      lambda_stat=(0.d0,0.d0)
	do k=1,nrr_div_c
	  k0=q_cube_in_sph_c(k)
	  if(k0==0) cycle
	  v=q_mesh_heg_c(:,k)-q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  v=-q_mesh_heg_c(:,k)+q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  lambda_stat(k)=b(k)+conjg(b(kq))
      enddo  !! over k
      deallocate(b)
      end