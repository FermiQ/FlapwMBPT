      subroutine k0_heg(ispin,q,ind_nu,nom,k0_pw,g_om,g_om_x)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,nom,ispin
      real*8, intent(in) :: q(*)
      complex*16, intent(in) :: g_om(nstar_c,nom,2),
     &                          g_om_x(nstar_c,nom,2)
	complex*16, intent(out) :: k0_pw(nrr_div_c,0:n_tau)
      integer :: i_omega,i_nu,k,k0,kq,kq0,i_tau
      real*8 :: v(3),de
      complex*16 :: cc
      complex*16, allocatable :: tmp1(:,:)
	i_nu=me_t*ndim3_nu+ind_nu-1
      allocate(tmp1(nrr_div_c,nom))
c ------------------------------------------------------------------
      tmp1=(0.d0,0.d0)
	do k=1,nrr_div_c
	  k0=q_cube_in_sph_c(k)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  v=q_mesh_heg_c(:,k)-q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  kq0=i_kref_c(kq0)
	  do i_omega=1,nom
	    tmp1(k,i_omega)=-g_om(k0,i_omega,1)*g_om(kq0,i_omega,2)
     &	                +g_om_x(k0,i_omega,1)*g_om_x(kq0,i_omega,2)
        enddo
      enddo  !! over i_omega
      call from_omega_nu_to_tau_nu_ab(nrr_div_c,ind_nu,tmp1,k0_pw,2)
c --------- Add the asymptotic correction --------------------------
	do k=1,nrr_div_c
	  k0=q_cube_in_sph_c(k)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  v=q_mesh_heg_c(:,k)-q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  kq0=i_kref_c(kq0)
	  de=e_star_xc(k0,ispin)-e_star_xc(kq0,ispin)
	  if(abs(de)<1.d-8.and.i_nu==0) then
	    do i_tau=0,n_tau
	      cc=gx_k_tau_heg_c(i_tau,k0,ispin)
     &	      *(betta_t*gx_k_tau_heg_c(n_tau,k0,ispin)
     &	       +tau_mesh(i_tau))/2.d0
           k0_pw(k,i_tau)=k0_pw(k,i_tau)+cc
          enddo
	  else
	    do i_tau=0,n_tau
	      cc=gx_k_tau_heg_c(i_tau,kq0,ispin)
     &        /dcmplx(de,-w_nu(i_nu))
            k0_pw(k,i_tau)=k0_pw(k,i_tau)+cc
          enddo
        endif
      enddo  !! over k
      deallocate(tmp1)
c ----- K0(k;tau) ----->  K0(r;tau) --------------------------
      do i_tau=0,n_tau
        call from_q_to_rr_heg_cube_red(k0_pw(1,i_tau))
      enddo  !! over i_tau
      end