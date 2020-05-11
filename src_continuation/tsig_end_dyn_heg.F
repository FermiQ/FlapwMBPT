      subroutine tsig_end_dyn_heg(nomnu,ind_omega,k1,t_pw,g_nu_om,
     &                            sig_om,wexa,wspl_geom,wspl_asy,wasy,
     &                            numtn)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: nomnu,ind_omega,k1,numtn
      real*8, intent(in) :: wexa(0:n_nu_exa,nstar_c),
     &                      wspl_geom(0:n_nu_geom+2,4,nstar_c),
     &                      wspl_asy(n_nu_asy+1,4,nstar_c),
     &                      wasy(nstar_c)
      complex*16, intent(in) :: g_nu_om(nstar_c,2,nomnu),
     &                          t_pw(nrr_div_c,numtn,2)
	complex*16, intent(inout) :: sig_om
      integer :: i_nu,kq,it,kq0,i,iq0,i1_omega,iq
      real*8 :: v(3),pib,w1
      complex*16 :: tt,tmp1
      real*8, allocatable :: w_nu_om(:,:,:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
	allocate(w_nu_om(nstar_c,2,nomnu))
	do i=1,nstar_c
	  do i_nu=1,nomnu
          w1=w_omega(i1_omega)-w_nu_adapt_omega(i_nu,i1_omega)-pib
          call interp_w_nu_spl(wexa(0,i),wspl_geom(0,1,i),
     &                         wspl_asy(1,1,i),wasy(i),
     &                         w_nu_om(i,1,i_nu),w1,1)
          w1=w_nu_adapt_omega(i_nu,i1_omega)
          call interp_w_nu_spl(wexa(0,i),wspl_geom(0,1,i),
     &                         wspl_asy(1,1,i),wasy(i),
     &                         w_nu_om(i,2,i_nu),w1,1)
        enddo
      enddo
c ------------------------------------------------------------------
	do iq=1,nrr_div_c
	  iq0=q_cube_in_sph_c(iq)
	  if(iq0==0) cycle
	  iq0=i_kref_c(iq0)
	  v=q_irr_c(:,k1)-q_mesh_heg_c(:,iq)
	  call zone1_number(v,rb0_c,nr_div_c,kq)
	  kq0=q_cube_in_sph_c(kq)
	  if(kq0==0) cycle
	  kq0=i_kref_c(kq0)
c -------------------------------------------------------------------
        tt=(0.d0,0.d0)     
	  do it=1,2
	    do i_nu=1,nomnu
            tmp1=w_nu_om(iq0,it,i_nu)*g_nu_om(kq0,it,i_nu)
     &                               *t_pw(iq,i_nu,it)
            tt=tt+tmp1*tau_om_from_nu_om_1(i_nu,0,ind_omega)
          enddo  !! over i_nu
        enddo
        sig_om=sig_om+tt/nqdiv_c/amega
      enddo
      deallocate(w_nu_om)
      end