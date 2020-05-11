      subroutine tsig0_end_dyn(nomnu,ind_omega,k1,t_sing,g_nu_om,sig_om,
     &                         ispin,nk,numtn)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nomnu,ind_omega,k1,ispin,nk,numtn
      complex*16, intent(in) :: g_nu_om(nbndf_bnd,nbndf_bnd,2,nomnu,
     &                                  npnt_c),
     &                          t_sing(nbndf_bnd,numtn,2)
      complex*16, intent(inout) :: sig_om(nbndf_bnd)
      integer :: i_nu,it,ii,i,j,i1_omega
      real*8 :: pib,w1,wasy,w_nu_om
      real*8, allocatable :: wexa(:),wspl_geom(:,:),wspl_asy(:,:)
      complex*16, allocatable :: tmp1(:),tt(:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
      allocate(tmp1(nk))
      allocate(wexa(0:n_nu_exa))
      allocate(wspl_geom(0:n_nu_geom+2,4))
      allocate(wspl_asy(n_nu_asy+1,4))
      call w_spl_loc(1,1,1,ws_head_nu,wexa,wspl_geom,wspl_asy,wasy)
c -------------------------------------------------------------------
      allocate(tt(nk))
      tt=(0.d0,0.d0)
      do it=1,2
        do i_nu=1,nomnu
          if(it==1) w1=w_omega(i1_omega)
     &                -w_nu_adapt_omega(i_nu,i1_omega)-pib
          if(it==2) w1=w_nu_adapt_omega(i_nu,i1_omega)
          call interp_w_nu_spl(wexa,wspl_geom,wspl_asy,wasy,w_nu_om,w1,
     &                         1)
          tmp1=(0.d0,0.d0)
          do j=1,nk
            do i=1,nk
              tmp1(i)=tmp1(i)+g_nu_om(i,j,it,i_nu,k1)*t_sing(j,i_nu,it)
            enddo
          enddo
          tmp1=tmp1*w_nu_om*q2aver_c
          tt=tt+tmp1*tau_om_from_nu_om_1(i_nu,0,ind_omega)
        enddo  !! over i_nu
      enddo
      deallocate(tmp1,wexa,wspl_geom,wspl_asy)
      w1=1.d0/sqrt(amega)/nqdiv_c
      do ii=1,nk
        sig_om(ii)=sig_om(ii)+w1*tt(ii)
      enddo
      deallocate(tt)
      end
