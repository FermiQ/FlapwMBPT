      subroutine t_heg_uniform(q_tau,t_tau)
	use atom_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: q_tau(nrr_div_c,ndim3_tau)
      real*8, intent(out) :: t_tau(nrr_div_c,ndim3_tau)
      integer :: i_nu1,iq1,k0,ind_tau,ind_nu
      real*8 :: w1
      real*8, allocatable :: t_nu(:,:),q_nu(:,:)
      complex*16, allocatable :: tm(:)
c -------- Transform Q(q';tau)[q;;nu] ---- > Q(q';nu')[q;;nu] ---------
      allocate(q_nu(nrr_div_c,ndim3_nu))
      call tau_to_nu_aa_g(q_nu,q_tau,nrr_div_c)
c ---------- Get T(q';nu')[q;;nu]   -------------------------------
      allocate(t_nu(nrr_div_c,ndim3_nu))
      t_nu=0.d0
	do iq1=1,nrr_div_c
	  k0=q_cube_in_sph_c(iq1)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  do ind_nu=1,ndim3_nu
	    i_nu1=me_t*ndim3_nu+ind_nu-1
          w1=w_q_nu_heg_c(i_nu1,k0)+v_q_heg_c(k0)
          t_nu(iq1,ind_nu)=w1*q_nu(iq1,ind_nu)*w1
c --- Subtract the static part ----------------------------------------
          t_nu(iq1,ind_nu)=t_nu(iq1,ind_nu)
     &                    -v_q_heg_c(k0)*q_nu(iq1,ind_nu)*v_q_heg_c(k0)
        enddo  !! over ind_nu
      enddo  !! over iq1
      deallocate(q_nu)
c ---------- Get T(q';tau)   -------------------------------
      call nu_to_tau_aa_g(t_nu,t_tau,nrr_div_c)
	deallocate(t_nu)
c --- Add static part ----------------------------------------
	do iq1=1,nrr_div_c
	  k0=q_cube_in_sph_c(iq1)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
        do ind_tau=1,ndim3_tau
          t_tau(iq1,ind_tau)=t_tau(iq1,ind_tau)
     &           +v_q_heg_c(k0)*q_tau(iq1,ind_tau)*v_q_heg_c(k0)
        enddo  !! over ind_tau
      enddo  !! over iq1
c ----- T(q';tau) ----->  T(r;tau) --------------------------
      allocate(tm(nrr_div_c))
      do ind_tau=1,ndim3_tau
	  tm=t_tau(:,ind_tau)
        call from_q_to_rr_heg_cube_red(tm)
	  t_tau(:,ind_tau)=tm
      enddo  !! over ind_tau
      deallocate(tm)
      end