      subroutine k0_heg_stat(ispin,q,k0_pw)
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
      real*8, intent(in) :: q(*)
	  complex*16, intent(out) :: k0_pw(nrr_div_c,2,ndim3_tau)
      integer :: i_omega,k,k0,kq,kq0,i_tau,ind_tau,it,ind_omega
      real*8 :: v(3),de
      complex*16 :: cc,ccx
      complex*16, allocatable :: tmp1(:,:,:),tmp(:,:),a_tau(:,:),
     &                           b_tau(:,:)
      allocate(tmp1(nrr_div_c,0:n_omega,2))
      allocate(tmp(nrr_div_c,ndim3_omega))
c -------------------------------------------------------------------
      tmp=(0.d0,0.d0)
	  do k=1,nrr_div_c
	    k0=i_kref_c(k)
	    v=q_sph_c(:,k)-q(1:3)
	    call zone1_number(v,rb0_c,nr_div_c,kq)
	    kq0=i_kref_c(kq)
	    do ind_omega=1,ndim3_omega
	      i_omega=me_t*ndim3_omega+ind_omega-1
	      cc=-g_k_omega_heg_c(i_omega,k0,ispin)
     & 	     *g_k_omega_heg_c(i_omega,kq0,ispin)
	      ccx=-gx_k_omega_heg_c(i_omega,k0,ispin)
     &	      *gx_k_omega_heg_c(i_omega,kq0,ispin)
	      tmp(k,ind_omega)=cc-ccx
        enddo
      enddo  !! over k
c --------- Transform to A and B form ---------------------------------
c --------- A=tmp1(+) + tmp1(-) ---------------------------------------
c --------- B=i{-tmp1(+) + tmp1(-)} -----------------------------------
      tmp1=(0.d0,0.d0)
	  do ind_omega=1,ndim3_omega
	    i_omega=me_t*ndim3_omega+ind_omega-1
	    do k=1,nrr_div_c
	      k0=i_kref_c(k)
	      v=-q_sph_c(:,k)+q(1:3)
	      call zone1_number(v,rb0_c,nr_div_c,kq)
	      kq0=i_kref_c(kq)
	      tmp1(k,i_omega,1)=tmp(k,ind_omega)+conjg(tmp(kq,ind_omega))
	      tmp1(k,i_omega,2)=(0.d0,1.d0)*(-tmp(k,ind_omega)
     &	                   +conjg(tmp(kq,ind_omega)))
	    enddo
	  enddo
	  if(nproc_t/=1) call dgop(tmp1,4*nrr_div_c*(1+n_omega),'  +',
     &	                         comm_t)
c --------- Get A(tau) and B(tau) -------------------------------------
      allocate(a_tau(nrr_div_c,ndim3_tau))
      allocate(b_tau(nrr_div_c,ndim3_tau))
	  call dgemm('n','n',2*nrr_div_c,ndim3_tau,n_omega+1,1.d0,
     &           tmp1(1,0,1),2*nrr_div_c,
     &           tau_from_omega(0,me_t*ndim3_tau,1),n_omega+1,0.d0,
     &           a_tau,2*nrr_div_c)
	  call dgemm('n','n',2*nrr_div_c,ndim3_tau,n_omega+1,1.d0,
     &           tmp1(1,0,2),2*nrr_div_c,
     &           tau_from_omega(0,me_t*ndim3_tau,2),n_omega+1,0.d0,
     &           b_tau,2*nrr_div_c)
c -------- Get tau and b-tau components of K0_PW ----------------------
      do ind_tau=1,ndim3_tau
        k0_pw(:,1,ind_tau)=0.25d0*(a_tau(:,ind_tau)+b_tau(:,ind_tau))
        k0_pw(:,2,ind_tau)=0.25d0*(b_tau(:,ind_tau)-a_tau(:,ind_tau))
      enddo
      deallocate(tmp,a_tau,b_tau,tmp1)
c --------- Add the asymptotic correction --------------------------
	  do k=1,nrr_div_c
	    k0=i_kref_c(k)
	    v=q_sph_c(:,k)-q(1:3)
	    call zone1_number(v,rb0_c,nr_div_c,kq)
	    kq0=i_kref_c(kq)
	    de=e_star_xc(kq0,ispin)-e_star_xc(k0,ispin)
	    do ind_tau=1,ndim3_tau
	      i_tau=me_t*ndim3_tau+ind_tau-1
	      do it=1,2
	        if(it==2) i_tau=n_tau-i_tau
	        if(abs(de)<1.d-8) then
	          cc=gx_k_tau_heg_c(i_tau,k0,ispin)
     &	        *(betta_t*gx_k_tau_heg_c(n_tau,k0,ispin)
     &	         +tau_mesh(i_tau))
	        else
	          cc=(gx_k_tau_heg_c(i_tau,k0,ispin)
     &	        -gx_k_tau_heg_c(i_tau,kq0,ispin))/de
	        endif
            k0_pw(k,it,ind_tau)=k0_pw(k,it,ind_tau)+cc
          enddo
        enddo
      enddo  !! over k
c ----- K(k;tau) ----->  K(r;tau) --------------------------
	  do ind_tau=1,ndim3_tau
	    do it=1,2
          call from_q_to_rr_heg_cube_red(k0_pw(1,it,ind_tau))
        enddo
      enddo  !! over ind_tau
      end
