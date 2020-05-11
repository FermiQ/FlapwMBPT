      subroutine t_heg(q,q_tau,t_tau,ind_nu,w2)
	use atom_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu
      real*8, intent(in) :: w2(nstar_c,2,n_nu2_max),q(*)
      complex*16, intent(in) :: q_tau(nrr_div_c,0:n_tau/2)
      complex*16, intent(out) :: t_tau(nrr_div_c,0:n_tau/2)
      integer :: i_tau,i_nu,i_nu1,nom_nu,iq1,kq,k0,kq0
      real*8 :: v(3)
      real*8, allocatable :: t_nu(:,:),bq(:,:,:),q_nu(:,:)
	i_nu=me_t*ndim3_nu+ind_nu-1
	nom_nu=num_nu2_adapt(i_nu)
      allocate(bq(nrr_div_c,0:n_tau/2,4))
c -------- Transform Q(q';tau)[q;;nu] ---- > Q(q';nu')[q;;nu] ---------
      do iq1=1,nrr_div_c
	  v=-q_mesh_heg_c(:,iq1)+q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
        do i_tau=0,n_tau/2
          bq(iq1,i_tau,1)=2.d0*real(q_tau(kq,i_tau))
          bq(iq1,i_tau,2)=2.d0*imag(q_tau(kq,i_tau))
          bq(iq1,i_tau,3)=2.d0*real(q_tau(iq1,i_tau))
          bq(iq1,i_tau,4)=-2.d0*imag(q_tau(iq1,i_tau))
        enddo
      enddo
      allocate(q_nu(nrr_div_c,nom_nu))
	call from_tau_nu_to_nu_nu_spl(nrr_div_c,ind_nu,nom_nu,q_nu,bq)
      deallocate(bq)
c ---------- Get T(q';nu')[q;;nu]   -------------------------------
      allocate(t_nu(nrr_div_c,nom_nu))
      t_nu=0.d0
	do iq1=1,nrr_div_c
	  k0=q_cube_in_sph_c(iq1)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  v=q_mesh_heg_c(:,iq1)-q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  kq0=i_kref_c(kq0)
        do i_nu1=1,nom_nu
          t_nu(iq1,i_nu1)=w2(k0,1,i_nu1)*q_nu(iq1,i_nu1)*w2(kq0,2,i_nu1)
        enddo  !! over i_nu1
      enddo  !! over iq1
      deallocate(q_nu)
c ---------- Get T(q';tau)[q;;nu]   -------------------------------
      call from_nu_nu_to_tau_nu_re(nrr_div_c,ind_nu,nom_nu,t_nu,t_tau)
      deallocate(t_nu)
c --- Add static part ----------------------------------------
	do iq1=1,nrr_div_c
	  k0=q_cube_in_sph_c(iq1)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  v=q_mesh_heg_c(:,iq1)-q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  kq0=i_kref_c(kq0)
        do i_tau=0,n_tau/2
        enddo  !! over i_tau
      enddo  !! over iq1
c ----- T(q';tau) ----->  T(r;tau) --------------------------
      do i_tau=0,n_tau/2
        call from_q_to_rr_heg_cube_red(t_tau(1,i_tau))
      enddo  !! over i_tau
      end