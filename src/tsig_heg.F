      subroutine tsig_heg(ispin,k_pw,t_pw,numtn,ind_omega,nomnu)
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ispin,numtn,ind_omega,nomnu
      complex*16, intent(in) :: k_pw(nrr_div_c,0:n_tau,2)
      complex*16, intent(out) :: t_pw(nrr_div_c,numtn,2)
      integer :: i_tau,i1_tau,ir,it,iq
      complex*16 :: cc
      complex*16, allocatable :: vxt(:,:),aq(:,:),t_nu(:)
      t_pw=(0.d0,0.d0)
      do it=1,2
        do i_tau=0,n_tau
          i1_tau=n_tau-i_tau
          do ir=1,nrr_div_c
            t_pw(ir,i_tau+1,it)=-g_r_tau_heg_c(ir,i1_tau,ispin)
     &                        *k_pw(ir,i_tau,it)
          enddo
        enddo    !! over i_tau
c ----- T(q;tau)[k;omega] ----------------------------
        do i_tau=1,n_tau+1
          call from_rr_to_q_heg_cube_red(t_pw(1,i_tau,it))
        enddo  !! over i_tau
      enddo   !! over it
c -------- Transform Vrt(tau;omega) ---- > Vrt(nu;omega) ------------
	allocate(vxt(0:n_tau,2))
      allocate(aq(0:n_tau/2,4))
      allocate(t_nu(nomnu))
	do iq=1,nrr_div_c
        do it=1,2
	    do i_tau=0,n_tau
            vxt(i_tau,it)=t_pw(iq,i_tau+1,it)
          enddo
        enddo
        do it=1,2
          cc=(-1)**it*(0.d0,1.d0)
          aq=(0.d0,0.d0)
          do i_tau=0,n_tau/2
            i1_tau=n_tau-i_tau
            aq(i_tau,1)=vxt(i_tau,3-it)+vxt(i1_tau,3-it)
            aq(i_tau,2)=vxt(i_tau,3-it)-vxt(i1_tau,3-it)
            aq(i_tau,2)=cc*aq(i_tau,2)
            aq(i_tau,3)=vxt(i_tau,it)+vxt(i1_tau,it)
            aq(i_tau,4)=vxt(i_tau,it)-vxt(i1_tau,it)
            aq(i_tau,4)=cc*aq(i_tau,4)
          enddo
          call from_tau_om_to_nu_om(2,ind_omega,nomnu,t_nu,aq)
          t_pw(iq,1:nomnu,it)=t_nu
        enddo
      enddo   !! over iq
      deallocate(vxt,aq,t_nu)
      end