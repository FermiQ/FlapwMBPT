      subroutine sigma_vrt_dyn_heg(aw,bw,sig,ind_nu)
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
      complex*16, intent(in) :: aw(nq_sph_c,0:n_tau),
     &                          bw(nq_sph_c,0:n_tau)
	real*8, intent(inout) :: sig(nq_sph_c,0:n_tau)
      integer :: i_tau
      real*8, allocatable :: s(:,:)
c ------ Lambda_dyn and  W_full = V + W_dyn   ---------------------
      allocate(s(nq_sph_c,0:n_tau))
	s=0.d0
      do i_tau=0,n_tau
	  s(:,i_tau)=s(:,i_tau)
     &        -tau_from_nu_d(ind_nu,0,1)*bw(:,i_tau)
     &	    -tau_from_nu_d(ind_nu,i_tau,1)*real(aw(:,i_tau))
     &	    +tau_from_nu_d(ind_nu,i_tau,2)*imag(aw(:,i_tau))
      enddo
      s=s*2
      sig=sig+s
      deallocate(s)
      end
      
      
      
      subroutine sigma_g_lambda_heg_spl(ispin,lambda_stat_tau,
     &                                  lambda_stat_nu_d2,sig,
     &                                  w_q_nuc_heg)
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
      real*8, intent(in) :: w_q_nuc_heg(ndim3_nu,nstar_c)
      complex*16, intent(in) :: lambda_stat_nu_d2(nrr_div_c,nspin,
     &                                            ndim3_nu,
     &                                            ndim_k_red(me_k+1)),
     &                          lambda_stat_tau(nrr_div_c,2,ndim3_tau,
     &                                          ndim_k_red(me_k+1))
	real*8, intent(out) :: sig(nq_sph_c,2,ndim3_tau)
      integer :: i_tau,k,ind_nu,ind_tau,it,ind_q,iq,kq,kq0
      real*8 :: fc,v(3)
      real*8, allocatable :: lw_tau(:,:,:)
      complex*16, allocatable :: lambda_stat_nu(:,:),lw(:,:)
      sig=0.d0
c ------ Lambda_stat   W_stat = V ---------------------------------
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        fc=wgt_c(iq)/amega*v_q_heg_c(iq)
        do ind_tau=1,ndim3_tau
          do it=1,2
            i_tau=me_t*ndim3_tau+ind_tau-1
            if(it==2) i_tau=n_tau-i_tau
	      do k=1,nq_sph_c
	        v=q_sph_c(:,k)+q_irr_c(:,iq)
	        call zone1_number(v,rb0_c,nr_div_c,kq)
	        kq0=q_cube_in_sph_c(kq)
	        if(kq0==0) cycle
	        kq0=i_kref_c(kq0)
	        sig(k,it,ind_tau)=sig(k,it,ind_tau)
     &	                     -lambda_stat_tau(kq,3-it,ind_tau,ind_q)
     &	                     *g_k_tau_heg_c(i_tau,kq0,ispin)*fc
            enddo  !! over k
          enddo  !! over it
        enddo  !! over ind_tau0
      enddo   !! over ind_q
c ------ Lambda_stat (order>1)   W_dyn ---------------------------------
      allocate(lambda_stat_nu(nrr_div_c,ndim3_nu))
      allocate(lw(nrr_div_c,ndim3_nu))
      allocate(lw_tau(nrr_div_c,2,ndim3_tau))
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
	  do ind_nu=1,ndim3_nu
	    lambda_stat_nu(:,ind_nu)=
     &      lambda_stat_nu_d2(:,ispin,ind_nu,ind_q)
	  enddo
        do ind_nu=1,ndim3_nu
          lw(:,ind_nu)=lambda_stat_nu(:,ind_nu)*w_q_nuc_heg(ind_nu,iq)
        enddo
        call nu_to_tau_c_spl(lw,lw_tau,nrr_div_c)
        fc=wgt_c(iq)/amega
        do ind_tau=1,ndim3_tau
          do it=1,2
            i_tau=me_t*ndim3_tau+ind_tau-1
            if(it==2) i_tau=n_tau-i_tau
	      do k=1,nq_sph_c
	        v=q_sph_c(:,k)+q_irr_c(:,iq)
	        call zone1_number(v,rb0_c,nr_div_c,kq)
	        kq0=q_cube_in_sph_c(kq)
	        if(kq0==0) cycle
	        kq0=i_kref_c(kq0)
	        sig(k,it,ind_tau)=sig(k,it,ind_tau)
     &	                     -lw_tau(kq,3-it,ind_tau)
     &	                     *g_k_tau_heg_c(i_tau,kq0,ispin)*fc
            enddo  !! over k
          enddo  !! over it
        enddo  !! over ind_tau
      enddo  !! over ind_q
      deallocate(lambda_stat_nu,lw,lw_tau)
	if(nproc_k/=1) call dgop(sig,2*nq_sph_c*ndim3_tau,'  +',comm_k)
      end
      
      
      