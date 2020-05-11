      subroutine k0_heg_uniform(ispin,k0_pw)
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
	real*8, intent(out) :: k0_pw(nrr_div_c,2,ndim3_tau)
      integer :: i_omega,k,k0,i_tau,ind_tau,it,ind_omega
      real*8 :: dd
      complex*16 :: cc,ccx
      complex*16, allocatable :: tmp(:,:),tm(:)
      allocate(tmp(nrr_div_c,ndim3_omega))
c -------------------------------------------------------------------
      tmp=(0.d0,0.d0)
	do k=1,nrr_div_c
	  k0=q_cube_in_sph_c(k)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  do ind_omega=1,ndim3_omega
	    i_omega=me_t*ndim3_omega+ind_omega-1
	    cc=-g_k_omega_heg_c(i_omega,k0,ispin)
     &	     *g_k_omega_heg_c(i_omega,k0,ispin)
	    ccx=-gx_k_omega_heg_c(i_omega,k0,ispin)
     &	     *gx_k_omega_heg_c(i_omega,k0,ispin)
	    tmp(k,ind_omega)=cc-ccx
        enddo
      enddo  !! over k
      call omega_to_tau_all(nrr_div_c,tmp,k0_pw)
      deallocate(tmp)
c --------- Add the asymptotic correction --------------------------
	do k=1,nrr_div_c
	  k0=q_cube_in_sph_c(k)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  do ind_tau=1,ndim3_tau
	    i_tau=me_t*ndim3_tau+ind_tau-1
	    do it=1,2
	      if(it==2) i_tau=n_tau-i_tau
	      dd=gx_k_tau_heg_c(i_tau,k0,ispin)
     &	      *(betta_t*gx_k_tau_heg_c(n_tau,k0,ispin)
     &	       +tau_mesh(i_tau))
            k0_pw(k,it,ind_tau)=k0_pw(k,it,ind_tau)+dd
          enddo
        enddo
      enddo  !! over k
c ----- K(k;tau) ----->  K(r;tau) --------------------------
      allocate(tm(nrr_div_c))
	do ind_tau=1,ndim3_tau
	  do it=1,2
	    tm=k0_pw(:,it,ind_tau)
          call from_q_to_rr_heg_cube_red(tm)
	    k0_pw(:,it,ind_tau)=tm
        enddo
      enddo  !! over ind_tau
      deallocate(tm)
      end