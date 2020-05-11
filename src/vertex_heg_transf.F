      subroutine vertex_heg_transf(q,ind_nu,nom,lambda_dyn_tau,
     &                             lambda_dyn)
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: nom,ind_nu
      real*8, intent(in) :: q(*)
      complex*16, intent(in) :: lambda_dyn_tau(nrr_div_c,0:n_tau)
      complex*16, intent(out) :: lambda_dyn(nrr_div_c,nom)
      integer :: i_tau,i1_tau,k,kq,k0,kq0
      real*8 :: v(3)
      complex*16, allocatable :: bk(:,:,:),ck(:,:)
      allocate(ck(nrr_div_c,0:n_tau))
      allocate(bk(nrr_div_c,0:n_tau/2,4))
c -------- Transform Vrt(R;tau) ---- > Vrt(K;tau) ----------
      ck=lambda_dyn_tau
      do i_tau=0,n_tau
        call from_rr_to_q_heg_cube_red(ck(1,i_tau))
      enddo  !! over i_tau
c -------- Transform Vrt(K;tau) ---- > Vrt(K;omega) ----------
      bk=(0.d0,0.d0)
	do k=1,nrr_div_c
	  k0=q_cube_in_sph_c(k)
	  if(k0==0) cycle
	  v=q_mesh_heg_c(:,k)-q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  v=-q_mesh_heg_c(:,k)+q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  do i_tau=0,n_tau/2
	    i1_tau=n_tau-i_tau
	    bk(k,i_tau,1)=conjg(ck(kq,i_tau))-conjg(ck(kq,i1_tau))
	    bk(k,i_tau,2)=conjg(ck(kq,i_tau))+conjg(ck(kq,i1_tau))
	    bk(k,i_tau,3)=ck(k,i_tau)-ck(k,i1_tau)
	    bk(k,i_tau,4)=ck(k,i_tau)+ck(k,i1_tau)
        enddo
      enddo
	call from_tau_nu_to_omega_nu_spl(nrr_div_c,ind_nu,nom,
     &                                 lambda_dyn,bk,2)
      deallocate(bk,ck)
      end