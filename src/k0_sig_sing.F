      subroutine k0_sig_sing(ind_omega,k1,k_sing,g_om_om,s2,nom,ispin)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,k1,s2,nom,ispin
      complex*16, intent(in) :: g_om_om(nbndf_bnd,nbndf_bnd,2,nom,
     &                                  npnt_c)
      complex*16, intent(out) :: k_sing(nbndf_bnd,0:n_tau,2)
      integer :: i_omega,it,nk,i1_omega,i_tau,kf0
      real*8 :: w1,pib,wasy,w_om_om
      real*8, allocatable :: wexa(:),wspl_geom(:,:),wspl_asy(:,:)
      complex*16, allocatable :: tmp1(:,:),tmp2(:,:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
      allocate(wexa(0:n_nu_exa))
      allocate(wspl_geom(0:n_nu_geom+2,4))
      allocate(wspl_asy(n_nu_asy+1,4))
c ------------------------------------------------------------------
      k_sing=(0.d0,0.d0)
      kf0=k_a_from_c(k1)
      nk=n_low_bnd(kf0,ispin)
      call w_spl_loc(1,1,1,ws_head_nu,wexa,wspl_geom,wspl_asy,wasy)
      allocate(tmp1(nk,n_omega1_max))
      allocate(tmp2(nk,0:n_tau))
      do it=1,2
        tmp1=(0.d0,0.d0)
        do i_omega=1,nom
          if(it==1) w1=w_om_adapt_omega(i_omega,i1_omega)+pib
          if(it==2) w1=w_omega(i1_omega)
     &                -w_om_adapt_omega(i_omega,i1_omega)
          call interp_w_nu_spl(wexa,wspl_geom,wspl_asy,wasy,w_om_om,w1,
     &                         1)
          tmp1(:,i_omega)=w_om_om*q2aver_c/sqrt(amega)
     &                   *g_om_om(1:nk,s2,it,i_omega,k1)
          if(it==1) tmp1(:,i_omega)=conjg(tmp1(:,i_omega))
        enddo  !! over i_omega
        call zgemm('n','n',nk,n_tau+1,nom,(1.d0,0.d0),tmp1,nk,
     &             tau_om_from_omega_om(1,0,ind_omega),n_omega1_max,
     &             (0.d0,0.d0),tmp2,nk)
        if(it==1) tmp2=conjg(tmp2)
        do i_tau=0,n_tau
          k_sing(1:nk,i_tau,it)=tmp2(:,i_tau)
        enddo
      enddo
      deallocate(tmp1,tmp2,wexa,wspl_geom,wspl_asy)
      end
