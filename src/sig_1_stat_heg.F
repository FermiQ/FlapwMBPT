      subroutine sig_1_stat_heg(it,ind_tau,t_pw,k1,ispin,sig_stat)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,ind_tau,k1,ispin
      complex*16, intent(inout) :: sig_stat
      complex*16, intent(out) :: t_pw(nrr_div_c)
      integer :: i_tau,i1_tau,iq,ir,iq0,kq,kq0
      real*8 :: v(3)
      complex*16, allocatable :: k_pw(:)
      i_tau=me_t*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau
      i1_tau=n_tau-i_tau
      allocate(k_pw(nrr_div_c))
      call k_sig_stat_heg(ispin,i_tau,k1,k_pw)
      do ir=1,nrr_div_c
        t_pw(ir)=-g_r_tau_heg_c(ir,i1_tau,ispin)*k_pw(ir)
      enddo
      call from_rr_to_q_heg_cube_red(t_pw)
	  do iq=1,nrr_div_c
	    iq0=q_cube_in_sph_c(iq)
	    if(iq0==0) cycle
	    iq0=i_kref_c(iq0)
	    v=q_irr_c(:,k1)-q_mesh_heg_c(:,iq)
	    call zone1_number(v,rb0_c,nr_div_c,kq)
	    kq0=q_cube_in_sph_c(kq)
	    if(kq0==0) cycle
	    kq0=i_kref_c(kq0)
        sig_stat=sig_stat+g_k_tau_heg_c(i_tau,kq0,ispin)*v_q_heg_c(iq0)
     &                                                  *t_pw(iq)
      enddo   !! over k
      deallocate(k_pw)
      sig_stat=sig_stat/nqdiv_c/amega
      end
