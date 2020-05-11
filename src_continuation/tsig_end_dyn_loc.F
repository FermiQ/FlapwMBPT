      subroutine tsig_end_dyn_loc(nomnu,ind_omega,t_pw,g_nu_om,
     &                            sig_om,wexa,wspl_geom,wspl_asy,
     &                            wasy,numtn,fif,n,np)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nomnu,ind_omega,numtn,n,np
      real*8, intent(in) :: wexa(np,np,0:n_nu_exa),
     &                      wspl_geom(np,np,0:n_nu_geom+2,4),
     &                      wspl_asy(np,np,n_nu_asy+1,4),
     &                      wasy(np,np),
     &                      fif(nrel,n,n,np)
      complex*16, intent(in) :: g_nu_om(n,n,2,nomnu),
     &                          t_pw(np,n,numtn,2)
      complex*16, intent(inout) :: sig_om(n)
      integer :: i_nu,it,ii,i,j,i1_omega
      real*8 :: pib,w1
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tt(:,:),
     &                           w_nu_om(:,:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
      allocate(tmp(np,n))
      allocate(tmp1(np,n))
      allocate(w_nu_om(np,np))
c -------------------------------------------------------------------
      allocate(tt(np,n))
      tt=(0.d0,0.d0)
      do it=1,2
        do i_nu=1,nomnu
          if(it==1) w1=w_omega(i1_omega)
     &                -w_nu_adapt_omega(i_nu,i1_omega)-pib
          if(it==2) w1=w_nu_adapt_omega(i_nu,i1_omega)
          call interp_w_nu_re_spl(wexa,wspl_geom,wspl_asy,wasy,
     &                            w_nu_om,w1,np,np)
          call zgemm('n','n',np,n,np,(1.d0,0.d0),w_nu_om,np,
     &               t_pw(1,1,i_nu,it),np,(0.d0,0.d0),tmp,np)
          call zgemm('n','t',np,n,n,(1.d0,0.d0),tmp,np,
     &               g_nu_om(1,1,it,i_nu),n,(0.d0,0.d0),tmp1,np)
          tt=tt+tmp1*tau_om_from_nu_om_1(i_nu,0,ind_omega)
        enddo  !! over i_nu
      enddo
      deallocate(tmp,tmp1,w_nu_om)
      do ii=1,n
        do j=1,n
          do i=1,np
            if(irel/=2) cc=fif(1,ii,j,i)
            if(irel==2) cc=dcmplx(fif(1,ii,j,i),fif(2,ii,j,i))
            sig_om(ii)=sig_om(ii)+cc*tt(i,j)
          enddo
        enddo
      enddo
      deallocate(tt)
      end
