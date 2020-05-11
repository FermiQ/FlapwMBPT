      subroutine k_sig_heg(ind_omega,k1,k_pw,g_om_om,nom,wexa,
     &                     wspl_geom,wspl_asy,wasy)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,k1,nom
      real*8, intent(in) :: wexa(0:n_nu_exa,nstar_c),
     &                      wspl_geom(0:n_nu_geom+2,4,nstar_c),
     &                      wspl_asy(n_nu_asy+1,4,nstar_c),
     &                      wasy(nstar_c)
      complex*16, intent(in) :: g_om_om(nstar_c,2,nom)
	complex*16, intent(out) :: k_pw(nrr_div_c,0:n_tau,2)
      integer :: i_omega,k,kq,it,k0,kq0,i1_omega,i,i_tau
      real*8 :: v(3),w1,pib
      real*8, allocatable :: w_om_om(:,:,:)
      complex*16, allocatable :: tmp1(:),tmp2(:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
	allocate(w_om_om(nstar_c,2,nom))
	do i=1,nstar_c
	  do i_omega=1,nom
          w1=w_om_adapt_omega(i_omega,i1_omega)+pib
          call interp_w_nu_spl(wexa(0,i),wspl_geom(0,1,i),
     &                         wspl_asy(1,1,i),wasy(i),
     &                         w_om_om(i,1,i_omega),w1,1)
          w1=w_omega(i1_omega)-w_om_adapt_omega(i_omega,i1_omega)
          call interp_w_nu_spl(wexa(0,i),wspl_geom(0,1,i),
     &                         wspl_asy(1,1,i),wasy(i),
     &                         w_om_om(i,2,i_omega),w1,1)
        enddo
      enddo
c ------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
      allocate(tmp1(nom)) 
      allocate(tmp2(0:n_tau))
	do k=1,nrr_div_c
	  k0=q_cube_in_sph_c(k)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  v=q_irr_c(:,k1)-q_mesh_heg_c(:,k)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  kq0=i_kref_c(kq0)   
	  do it=1,2
	    tmp1=(0.d0,0.d0)
	    do i_omega=1,nom
	      tmp1(i_omega)=g_om_om(k0,it,i_omega)*w_om_om(kq0,it,i_omega)
          enddo  !! over i_omega
          if(it==1) tmp1=conjg(tmp1)
          call zgemm('n','n',1,n_tau+1,nom,(1.d0,0.d0),tmp1,1,
     &	           tau_om_from_omega_om(1,0,ind_omega),n_omega1_max,
     &               (0.d0,0.d0),tmp2,1)
          if(it==1) tmp2=conjg(tmp2)
          k_pw(k,:,it)=tmp2
        enddo
      enddo   !! over k
      deallocate(tmp1,tmp2,w_om_om)
c ----- K(k;tau) ----->  K(r;tau) --------------------------
      do it=1,2
        do i_tau=0,n_tau
          call from_q_to_rr_heg_cube_red(k_pw(1,i_tau,it))
        enddo  !! over i_tau
      enddo
      end
      
      
      
      
      subroutine k_sig_stat_heg(ispin,i_tau,k1,k_pw)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ispin,i_tau,k1
	complex*16, intent(out) :: k_pw(nrr_div_c)
      integer :: k,kq,k0,kq0
      real*8 :: v(3)
      k_pw=(0.d0,0.d0)
	do k=1,nrr_div_c
	  k0=q_cube_in_sph_c(k)
	  if(k0==0) cycle
	  k0=i_kref_c(k0)
	  v=q_irr_c(:,k1)-q_mesh_heg_c(:,k)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  kq0=i_kref_c(kq0)   
	  k_pw(k)=g_k_tau_heg_c(i_tau,k0,ispin)*v_q_heg_c(kq0)
      enddo   !! over k
c ----- K(k;tau) ----->  K(r;tau) --------------------------
      call from_q_to_rr_heg_cube_red(k_pw)
      end