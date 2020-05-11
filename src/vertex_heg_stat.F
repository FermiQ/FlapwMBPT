      subroutine vertex_heg_stat_spl(lambda_stat_tau,
     &                               lambda_stat_nu_d2,
     &                               lambda_stat_nu_d1)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      complex*16, intent(in) :: lambda_stat_nu_d2(nrr_div_c,nspin,
     &                                            ndim3_nu,
     &                                            ndim_k_red(me_k+1))
      complex*16, intent(out) :: lambda_stat_tau(nrr_div_c,2,ndim3_tau,
     &                                           ndim_k_red(me_k+1),
     &                                           nspin),
     &                           lambda_stat_nu_d1(nrr_div_c,ndim3_nu,
     &                                             ndim_k_red(me_k+1),
     &                                             nspin)
      integer :: ind_tau,it,i_tau,i,ind_q,iq,k,k0,kq,kq0,ind_nu,
     &           ispin
      real*8 :: v(3)
      real*8, allocatable :: d_lambda_tau(:,:,:)
      complex*16, allocatable :: d_lambda(:,:),k0t(:),d0_lambda(:,:)
      allocate(k0t(nrr_div_c))
      allocate(d0_lambda(nrr_div_c,0:n_nu))
      allocate(d_lambda(nrr_div_c,ndim3_nu))
      allocate(d_lambda_tau(nrr_div_c,2,ndim3_tau))
c -------- First Order Static Vertex (Tau) --------------------------
      lambda_stat_tau=(0.d0,0.d0)
      do ispin=1,nspin
        do ind_q=1,ndim_k_red(me_k+1)
          iq=n_mpi_k_red(me_k+1)+ind_q
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              k0t=(0.d0,0.d0)
	        do k=1,nrr_div_c
	          k0=q_cube_in_sph_c(k)
	          if(k0==0) cycle
	          k0=i_kref_c(k0)
	          v=q_mesh_heg_c(:,k)-q_irr_c(:,iq)
	          call zone1_number(v,rb0_c,nr_div_c,kq)
	          kq0=q_cube_in_sph_c(kq)
	          if(kq0==0) cycle
	          kq0=i_kref_c(kq0)
	          k0t(k)=g_k_tau_heg_c(i_tau,k0,ispin)
     &  	            *g_k_tau_heg_c(n_tau-i_tau,kq0,ispin)
              enddo
              call from_q_to_rr_heg_cube_red(k0t)
              do i=1,nrr_div_c
                k0t(i)=v_r_heg_c(i)*k0t(i)
              enddo
              call from_rr_to_q_heg_cube_red(k0t)
              lambda_stat_tau(:,it,ind_tau,ind_q,ispin)=k0t
            enddo
          enddo
          call tau_to_nu_gen_spl(lambda_stat_nu_d1(1,1,ind_q,ispin),
     &                           lambda_stat_tau(1,1,1,ind_q,ispin),
     &                           nrr_div_c)
c -------- Higher Order Static Vertex (Tau) --------------------------
	    do ind_nu=1,ndim3_nu
	      d_lambda(:,ind_nu)=lambda_stat_nu_d2(:,ispin,ind_nu,ind_q)
	    enddo
          call nu_to_tau_c_spl(d_lambda,d_lambda_tau,nrr_div_c)
          do ind_tau=1,ndim3_tau
            do it=1,2
              lambda_stat_tau(:,it,ind_tau,ind_q,ispin)=
     &        lambda_stat_tau(:,it,ind_tau,ind_q,ispin)
     &       +d_lambda_tau(:,it,ind_tau)
            enddo
          enddo
        enddo   !! over ind_q
      enddo
      deallocate(k0t,d_lambda,d_lambda_tau,d0_lambda)
      end