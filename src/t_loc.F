      subroutine t_loc(qt,vloc,ind_nu,nomnu,wloc,np)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,nomnu,np
      real*8, intent(in) :: wloc(n_pbmtm_red,n_pbmtm_red,0:n_nu),
     &                      vloc(n_pbmtm_red,n_pbmtm_red)
      complex*16, intent(inout) :: qt(np,np,2,2,0:n_tau/2)
      integer :: j,i,nn,i_nu1,in,i_tau,i_nu,i0_tau,it,il
      real*8 :: p2b,w1
      complex*16 :: cc
      real*8, allocatable :: wexa(:,:,:),wspl_geom(:,:,:,:),
     &                       wspl_asy(:,:,:,:),wasy(:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),
     &                           t_nu(:,:,:),aq(:,:,:,:),
     &                           q_pw(:,:,:,:,:),w_nu_nu(:,:),
     &                           q_tau(:,:,:)
      p2b=(pi+pi)/betta_t
      i_nu=me_t*ndim3_nu+ind_nu-1
      allocate(tmp(np,np))
      allocate(tmp1(np,np))
      allocate(q_pw(np,np,2,0:n_tau/2,2))
      allocate(w_nu_nu(np,np))
      allocate(wexa(np,np,0:n_nu_exa))
      allocate(wspl_geom(np,np,0:n_nu_geom+2,4))
      allocate(wspl_asy(np,np,n_nu_asy+1,4))
      allocate(wasy(np,np))
      call w_spl_loc(np,n_pbmtm_red,np,wloc,wexa,wspl_geom,wspl_asy,
     &               wasy)
      nn=np**2
      do i0_tau=0,n_tau/2
        do il=1,2
          do it=1,2
            q_pw(:,:,it,i0_tau,il)=qt_loc(:,:,il,it,i0_tau)
          enddo
        enddo
      enddo
c -------- Transform Q(tau;nu) ---- > Q(nu';nu) ---------
      do in=1,2
        cc=(-1)**in*(0.d0,1.d0)
        allocate(aq(np,np,0:n_tau/2,4))
        aq=(0.d0,0.d0)
        do i_tau=0,n_tau/2
          do j=1,np
            do i=1,np
              aq(i,j,i_tau,1)=q_pw(i,j,1,i_tau,3-in)
     &                       +q_pw(i,j,2,i_tau,3-in)
              aq(i,j,i_tau,2)=q_pw(i,j,1,i_tau,3-in)
     &                       -q_pw(i,j,2,i_tau,3-in)
              aq(i,j,i_tau,2)=cc*aq(i,j,i_tau,2)
              aq(i,j,i_tau,3)=q_pw(i,j,1,i_tau,in)
     &                       +q_pw(i,j,2,i_tau,in)
              aq(i,j,i_tau,4)=q_pw(i,j,1,i_tau,in)
     &                       -q_pw(i,j,2,i_tau,in)
              aq(i,j,i_tau,4)=cc*aq(i,j,i_tau,4)
            enddo
          enddo
        enddo
        allocate(t_nu(np,np,nomnu))
        call from_tau_nu_to_nu_nu_spl(2*nn,ind_nu,nomnu,t_nu,aq)
        deallocate(aq)
c ---------- Get T_q_nu and subtract static asymptotics ----------
        do i_nu1=1,nomnu
          if(in==1) w1=w_nu_adapt_nu(i_nu1,i_nu)-w_nu(i_nu)
          if(in==2) w1=w_nu_adapt_nu(i_nu1,i_nu)
          call interp_w_nu_re_spl(wexa,wspl_geom,wspl_asy,wasy,w_nu_nu,
     &                            w1,np,np)
          w_nu_nu(1:np,1:np)=w_nu_nu(1:np,1:np)+vloc(1:np,1:np)
          call zgemm('n','n',np,np,np,(1.d0,0.d0),w_nu_nu,np,
     &               t_nu(1,1,i_nu1),np,(0.d0,0.d0),tmp,np)
          if(in==1) w1=w_nu_adapt_nu(i_nu1,i_nu)
          if(in==2) w1=w_nu_adapt_nu(i_nu1,i_nu)-w_nu(i_nu)
          call interp_w_nu_re_spl(wexa,wspl_geom,wspl_asy,wasy,w_nu_nu,
     &                            w1,np,np)
          w_nu_nu(1:np,1:np)=w_nu_nu(1:np,1:np)+vloc(1:np,1:np)
          call zgemm('n','n',np,np,np,(1.d0,0.d0),tmp,np,w_nu_nu,np,
     &               (0.d0,0.d0),t_nu(1,1,i_nu1),np)
        enddo  !! over i_nu1
c ---------------------------------------------------------------
        if(in==1) t_nu=conjg(t_nu)
        allocate(q_tau(np,np,0:n_tau))
        call zgemm('n','n',nn,n_tau+1,nomnu,(1.d0,0.d0),t_nu,nn,
     &             tau_nu_from_nu_nu_1(1,0,ind_nu),n_nu2_max,
     &             (0.d0,0.d0),q_tau,nn)
        deallocate(t_nu)
        if(in==1) q_tau=conjg(q_tau)
        do i0_tau=0,n_tau/2
          do it=1,2
            i_tau=i0_tau
            if(it==2) i_tau=n_tau-i0_tau
            qt_loc(:,:,in,it,i0_tau)=q_tau(:,:,i_tau)
          enddo
        enddo
        deallocate(q_tau)
      enddo  !! over in
      deallocate(tmp,tmp1,q_pw,w_nu_nu,wexa,wspl_geom,wspl_asy,wasy)
      end
