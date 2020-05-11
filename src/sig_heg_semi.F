      subroutine sig_heg_semi(ispin,lambda_stat_nu_d1,sig)
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
      complex*16, intent(in) :: lambda_stat_nu_d1(nrr_div_c,ndim3_nu,
     &                                            ndim_k_red(me_k+1))
      real*8, intent(out) :: sig(nq_sph_c,2,ndim3_tau)
      integer :: ind_tau,it,ind_nu,k,i_tau,i_nu,ind_q,iq,kq,kq0
      real*8 :: fc,v(3)
      real*8, allocatable :: l_tau(:,:,:)
      complex*16, allocatable :: l_nu(:,:)
      allocate(l_nu(nrr_div_c,ndim3_nu))
      allocate(l_tau(nrr_div_c,2,ndim3_tau))
      sig=0.d0
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          l_nu(:,ind_nu)=lambda_stat_nu_d1(:,ind_nu,ind_q)
     &                  *w_q_nu_heg_c(i_nu,iq)
        enddo
        call nu_to_tau_c_spl(l_nu,l_tau,nrr_div_c)
c ------------------------------------------------------------------
        fc=2.d0*wgt_c(iq)/amega
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
	        sig(k,it,ind_tau)=sig(k,it,ind_tau)-l_tau(kq,3-it,ind_tau)
     &	                     *g_k_tau_heg_c(i_tau,kq0,ispin)*fc
            enddo  !! over k
          enddo  !! over it
        enddo  !! over ind_tau
      enddo   !! over ind_q
      deallocate(l_nu,l_tau)
	if(nproc_k/=1) call dgop(sig,2*nq_sph_c*ndim3_tau,'  +',comm_k)
      end