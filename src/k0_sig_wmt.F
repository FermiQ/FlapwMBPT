      subroutine k0_sig_wmt(ind_omega,jatom,k_pw,g_om_om,n,n1,np1,fif,
     &                      s2,nom,ispin,wloc)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,jatom,s2,nom,ispin,n,n1,np1
      real*8, intent(in) :: wloc(n_pbmtm_red,n_pbmtm_red,0:n_nu),
     &                      fif(nrel,maxel_red,maxel_red,n_pbmtm_red)
      complex*16, intent(in) :: g_om_om(n,n1,2,nom)
      complex*16, intent(out) :: k_pw(n,np1,0:n_tau,2)
      integer :: i_omega,nn,it,i1_omega,i_tau,i2,jsort,l,l2
      real*8 :: w1,pib
      complex*16 :: cc
      real*8, allocatable :: wexa(:,:,:),wspl_geom(:,:,:,:),
     &                       wspl_asy(:,:,:,:),wasy(:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:),tmp2(:,:,:),
     &                           w_om_om(:,:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
      allocate(tmp(n,np1))
      allocate(w_om_om(n_pbmtm_red,n_pbmtm_red))
      allocate(wexa(n_pbmtm_red,n_pbmtm_red,0:n_nu_exa))
      allocate(wspl_geom(n_pbmtm_red,n_pbmtm_red,0:n_nu_geom+2,4))
      allocate(wspl_asy(n_pbmtm_red,n_pbmtm_red,n_nu_asy+1,4))
      allocate(wasy(n_pbmtm_red,n_pbmtm_red))
c ------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
      nn=n*np1
      jsort=is(jatom)
      call w_spl_loc_red(jatom,wloc,wexa,wspl_geom,wspl_asy,wasy)
      allocate(tmp1(n,np1,n_omega1_max))
      allocate(tmp2(n,np1,0:n_tau))
      do it=1,2
        tmp1=(0.d0,0.d0)
        do i_omega=1,nom
          if(it==1) w1=w_om_adapt_omega(i_omega,i1_omega)+pib
          if(it==2) w1=w_omega(i1_omega)
     &                -w_om_adapt_omega(i_omega,i1_omega)
          call interp_w_nu_re_spl(wexa,wspl_geom,wspl_asy,wasy,
     &                            w_om_om,w1,n_pbmt0_red(jsort),
     &                            n_pbmtm_red)
          tmp=(0.d0,0.d0)
          do i2=1,np1
            do l2=1,n1
              if(irel/=2) cc=fif(1,l2,s2,i2)
              if(irel==2) cc=dcmplx(fif(1,l2,s2,i2),fif(2,l2,s2,i2))
              do l=1,n
                tmp(l,i2)=tmp(l,i2)+g_om_om(l,l2,it,i_omega)*cc
              enddo
            enddo
          enddo
          call zgemm('n','t',n,np1,np1,(1.d0,0.d0),tmp,n,
     &               w_om_om,n_pbmtm_red,(0.d0,0.d0),
     &               tmp1(1,1,i_omega),n)
          if(it==1) tmp1(:,:,i_omega)=conjg(tmp1(:,:,i_omega))
        enddo  !! over i_omega
        call zgemm('n','n',nn,n_tau+1,nom,(1.d0,0.d0),tmp1,nn,
     &             tau_om_from_omega_om(1,0,ind_omega),n_omega1_max,
     &             (0.d0,0.d0),tmp2,nn)
        if(it==1) tmp2=conjg(tmp2)
        do i_tau=0,n_tau
          k_pw(1:n,1:np1,i_tau,it)=tmp2(:,:,i_tau)
        enddo
      enddo  !! over it
      deallocate(tmp1,tmp2,tmp,w_om_om,wexa,wspl_geom,wspl_asy,wasy)
      end
