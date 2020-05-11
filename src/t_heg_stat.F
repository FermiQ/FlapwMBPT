      subroutine t_heg_stat(q,q_tau,t_tau)
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
      complex*16, intent(in) :: q_tau(nrr_div_c,ndim3_tau)
      complex*16, intent(out) :: t_tau(nrr_div_c,ndim3_tau)
      integer :: i_nu1,iq1,kq,k0,kq0,ind_tau,k,ind_nu
      real*8 :: v(3),w1,w2
      real*8, allocatable :: t_nu(:,:),bq(:,:,:),q_nu(:,:),tmp(:,:),
     &                       tmp1(:,:,:),a_tau(:,:),b_tau(:,:)
      allocate(bq(nrr_div_c,ndim3_tau,2))
c -------- Transform Q(q';tau)[q;;nu] ---- > Q(q';nu')[q;;nu] ---------
      do ind_tau=1,ndim3_tau
        bq(:,ind_tau,1)=real(q_tau(:,ind_tau))
        bq(:,ind_tau,2)=imag(q_tau(:,ind_tau))
      enddo
      bq=2.d0*bq
      allocate(q_nu(nrr_div_c,0:n_nu))
      allocate(tmp(nrr_div_c,0:n_nu))
c ---------------------- COS(v')t contribution ------------------------
	call dgemm('n','n',nrr_div_c,n_nu+1,ndim3_tau,1.d0,bq(1,1,1),
     &	         nrr_div_c,ab_nu_from_tau(1,0,1),ndim3_tau,0.d0,tmp,
     &           nrr_div_c)
c ---------------------- SIN(v')t contribution ------------------------
	call dgemm('n','n',nrr_div_c,n_nu+1,ndim3_tau,-1.d0,bq(1,1,2),
     &	         nrr_div_c,ab_nu_from_tau(1,0,2),ndim3_tau,1.d0,tmp,
     &           nrr_div_c)
      deallocate(bq)
      deallocate(tmp)
c ---------- Get T(q';nu')[q;;nu]   -------------------------------
      allocate(t_nu(nrr_div_c,ndim3_nu))
      t_nu=0.d0
	do iq1=1,nrr_div_c
	  k0=i_kref_c(iq1)
	  v=q_sph_c(:,iq1)-q(1:3)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=i_kref_c(kq)
	  do ind_nu=1,ndim3_nu
	    i_nu1=me_t*ndim3_nu+ind_nu-1
          w1=w_q_nu_heg(i_nu1,k0)+v_q_heg(k0)
          w2=w_q_nu_heg(i_nu1,kq0)+v_q_heg(kq0)
          t_nu(iq1,ind_nu)=w1*q_nu(iq1,i_nu1)*w2
        enddo  !! over i_nu1
      enddo  !! over iq1
      deallocate(q_nu)
c --------- Transform to A and B form ---------------------------------
c --------- A=tmp1(+) + tmp1(-) ---------------------------------------
c --------- B=-tmp1(+) + tmp1(-) -----------------------------------
      allocate(tmp1(nrr_div_c,0:n_nu,2))
      tmp1=0.d0
	do ind_nu=1,ndim3_nu
	  i_nu1=me_t*ndim3_nu+ind_nu-1
	  do k=1,nrr_div_c
	    k0=i_kref_c(k)
	    v=-q_sph_c(:,k)+q(1:3)
	    call zone1_number(v,rb0_c,nr_div_c,kq)
	    kq0=i_kref_c(kq)
	    tmp1(k,i_nu1,1)=t_nu(k,ind_nu)+t_nu(kq,ind_nu)
	    tmp1(k,i_nu1,2)=-t_nu(k,ind_nu)+t_nu(kq,ind_nu)
	  enddo
	enddo
	deallocate(t_nu)
	if(nproc_t/=1) call dgop(tmp1,2*nrr_div_c*(1+n_nu),'  +',
     &	                         comm_t)
	allocate(a_tau(nrr_div_c,ndim3_tau))
	allocate(b_tau(nrr_div_c,ndim3_tau))
c ---------- Get T(q';tau)[q]   -------------------------------
	call dgemm('n','n',nrr_div_c,ndim3_tau,n_nu+1,1.d0,tmp1(1,0,1),
     &	   nrr_div_c,tau_from_nu(0,1,1),n_nu+1,0.d0,a_tau,nrr_div_c)
	call dgemm('n','n',nrr_div_c,ndim3_tau,n_nu+1,1.d0,tmp1(1,0,2),
     &	nrr_div_c,tau_from_nu(0,1,2),n_nu+1,0.d0,b_tau,nrr_div_c)
      t_tau=(a_tau+(0.d0,1.d0)*b_tau)/4
      deallocate(a_tau,b_tau)
c ----- T(q';tau) ----->  T(r;tau) --------------------------
      do ind_tau=1,ndim3_tau
        call from_q_to_rr_heg_cube_red(t_tau(1,ind_tau))
      enddo  !! over ind_tau
      end