      subroutine sig_1_semi_heg(t_pw,k1,ispin,sig_semi)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k1,ispin
      complex*16, intent(in) :: t_pw(nrr_div_c,2,ndim3_tau)
      real*8, intent(inout) :: sig_semi(2,ndim3_tau)
      integer :: iq,iq0,kq,kq0,it,i_nu,i_tau,ind_nu,ind_tau
      real*8 :: v(3)
      complex*16, allocatable :: t_nu(:,:,:),wtt(:,:,:)
      allocate(t_nu(nrr_div_c,ndim3_nu,2))
      call tau_to_nu_full(t_nu,t_pw,nrr_div_c)
      do iq=1,nrr_div_c
	  iq0=q_cube_in_sph_c(iq)
	  if(iq0==0) then
          t_nu(iq,:,:)=0.d0
          cycle
        endif
	  iq0=i_kref_c(iq0)
        do it=1,2
          do ind_nu=1,ndim3_nu
            i_nu=me_t*ndim3_nu+ind_nu-1
            t_nu(iq,ind_nu,it)=t_nu(iq,ind_nu,it)*w_q_nu_heg_c(i_nu,iq0)
          enddo
        enddo
      enddo
      allocate(wtt(nrr_div_c,2,ndim3_tau))
      call nu_to_tau_cc(t_nu,wtt,nrr_div_c)
      deallocate(t_nu)
	do iq=1,nrr_div_c
	  iq0=q_cube_in_sph_c(iq)
	  if(iq0==0) cycle
	  iq0=i_kref_c(iq0)
	  v=q_irr_c(:,k1)-q_mesh_heg_c(:,iq)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  kq0=i_kref_c(kq0)
        do ind_tau=1,ndim3_tau
          i_tau=me_t*ndim3_tau+ind_tau-1
          do it=1,2
            if(it==2) i_tau=n_tau-i_tau
	      sig_semi(it,ind_tau)=sig_semi(it,ind_tau)
     &                          +g_k_tau_heg_c(i_tau,kq0,ispin)
     &                          *wtt(iq,it,ind_tau)
          enddo
        enddo
      enddo   !! over k
      deallocate(wtt)
      sig_semi=sig_semi/nqdiv_c/amega
      end