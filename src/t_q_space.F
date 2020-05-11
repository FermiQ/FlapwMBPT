      subroutine t_q_space(qt_pw,v2_q,ind_nu,iq,nomnu,wnu)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,iq,nomnu
      real*8, intent(in) :: wnu(n_pbtot_red,n_pbtot_red,0:n_nu,npnt_c)
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c)
      complex*16, intent(inout) :: qt_pw(n_pbtot_red,n_pbtot_red,
     &                                   nqdiv_c,2,2,0:n_tau/2)
      integer :: j,i,nn,iq1,iq2,i_nu1,iq3,iq10,np,iq30,np3,iq20,np2,in,
     &           i_tau,ka,ka3,ka2,i_nu,i0_tau,it,il
      real*8 :: v(3),p2b,w1
      complex*16 :: cc
      real*8, allocatable :: wexa(:,:,:,:),wspl_geom(:,:,:,:,:),
     &                       wspl_asy(:,:,:,:,:),wasy(:,:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),
     &                           t_nu(:,:,:),aq(:,:,:,:),
     &                           q_pw(:,:,:,:,:),w_nu_nu(:,:),
     &                           q_tau(:,:,:)
	  p2b=(pi+pi)/betta_t
      i_nu=me_t*ndim3_nu+ind_nu-1
      allocate(tmp(n_pbtot_red,n_pbtot_red))
      allocate(tmp1(n_pbtot_red,n_pbtot_red))
      allocate(q_pw(n_pbtot_red,n_pbtot_red,2,0:n_tau/2,2))
      allocate(w_nu_nu(n_pbtot_red,n_pbtot_red))
      allocate(wexa(n_pbtot_red,n_pbtot_red,0:n_nu_exa,2))
      allocate(wspl_geom(n_pbtot_red,n_pbtot_red,0:n_nu_geom+2,4,2))
      allocate(wspl_asy(n_pbtot_red,n_pbtot_red,n_nu_asy+1,4,2))
      allocate(wasy(n_pbtot_red,n_pbtot_red,2))
      nn=n_pbtot_red**2
      do iq1=1,nqdiv_c
	    iq10=i_kref_c(iq1)
        ka=k_a_from_c(iq10)
	    np=n_pbmt_red+nplwgw_red(ka)
        v=pnt_c(:,iq)-pnt_c(:,iq1)
        call zone1_number(v,rb0_c,ndiv_c,iq3)  ! iq3 - index in BZ_1
        iq3=index_k1_c(iq3)  ! Now iq3 - index in BZ_0
	    iq30=i_kref_c(iq3)
        ka3=k_a_from_c(iq30)
	    np3=n_pbmt_red+nplwgw_red(ka3)
        v=pnt_c(:,iq1)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,iq2)  ! iq2 - index in BZ_1
        iq2=index_k1_c(iq2)  ! Now iq2 - index in BZ_0
	    iq20=i_kref_c(iq2)
        ka2=k_a_from_c(iq20)
	    np2=n_pbmt_red+nplwgw_red(ka2)
        do i0_tau=0,n_tau/2
          do il=1,2
            do it=1,2
	          q_pw(:,:,it,i0_tau,il)=qt_pw(:,:,iq1,il,it,i0_tau)
            enddo
          enddo
        enddo
        call w_spl_q_red(iq1,wnu,wexa(1,1,0,1),wspl_geom(1,1,0,1,1),
     &                   wspl_asy(1,1,1,1,1),wasy(1,1,1))
        call w_spl_q_red(iq2,wnu,wexa(1,1,0,2),wspl_geom(1,1,0,1,2),
     &                   wspl_asy(1,1,1,1,2),wasy(1,1,2))
c -------- Transform Q(tau;nu) ---- > Q(nu';nu) ---------
        do in=1,2
          cc=(-1)**in*(0.d0,1.d0)
          allocate(aq(np,np3,0:n_tau/2,4))
          aq=(0.d0,0.d0)
          do i_tau=0,n_tau/2
            do j=1,np3
              do i=1,np
                aq(i,j,i_tau,1)=q_pw(i,j,1,i_tau,3-in)
     &                         +q_pw(i,j,2,i_tau,3-in)
                aq(i,j,i_tau,2)=q_pw(i,j,1,i_tau,3-in)
     &                         -q_pw(i,j,2,i_tau,3-in)
                aq(i,j,i_tau,2)=cc*aq(i,j,i_tau,2)
                aq(i,j,i_tau,3)=q_pw(i,j,1,i_tau,in)
     &                         +q_pw(i,j,2,i_tau,in)
                aq(i,j,i_tau,4)=q_pw(i,j,1,i_tau,in)
     &                         -q_pw(i,j,2,i_tau,in)
                aq(i,j,i_tau,4)=cc*aq(i,j,i_tau,4)
              enddo
            enddo
          enddo
          allocate(t_nu(np,np3,nomnu))
          call from_tau_nu_to_nu_nu_spl(2*np*np3,ind_nu,nomnu,t_nu,aq)
          deallocate(aq)
c ---------- Get T_q_nu and subtract static asymptotics ----------
          do i_nu1=1,nomnu
            if(in==1) w1=w_nu_adapt_nu(i_nu1,i_nu)-w_nu(i_nu)
            if(in==2) w1=w_nu_adapt_nu(i_nu1,i_nu)
            call interp_w_nu_sq_spl(wexa(1,1,0,1),wspl_geom(1,1,0,1,1),
     &                              wspl_asy(1,1,1,1,1),wasy(1,1,1),
     &                              w_nu_nu,w1,n_pbtot_red)
            w_nu_nu=w_nu_nu+v2_q(:,:,iq1)
	        call zgemm('n','n',np,np2,np,(1.d0,0.d0),w_nu_nu,
     &	               n_pbtot_red,t_nu(1,1,i_nu1),np,(0.d0,0.d0),tmp,
     &                 n_pbtot_red)
c	        call zgemm('n','n',np,np2,np,(1.d0,0.d0),v2_q(1,1,iq1),
c     &	               n_pbtot_red,t_nu(1,1,i_nu1),np,(0.d0,0.d0),
c     &                 tmp1,n_pbtot_red)
            if(in==1) w1=w_nu_adapt_nu(i_nu1,i_nu)
            if(in==2) w1=w_nu_adapt_nu(i_nu1,i_nu)-w_nu(i_nu)
            call interp_w_nu_sq_spl(wexa(1,1,0,2),wspl_geom(1,1,0,1,2),
     &                              wspl_asy(1,1,1,1,2),wasy(1,1,2),
     &                              w_nu_nu,w1,n_pbtot_red)
            w_nu_nu=w_nu_nu+v2_q(:,:,iq2)
	        call zgemm('n','n',np,np2,np2,(1.d0,0.d0),tmp,n_pbtot_red,
     &	               w_nu_nu,n_pbtot_red,(0.d0,0.d0),
     &                 t_nu(1,1,i_nu1),np)
c	        call zgemm('n','n',np,np2,np2,(-1.d0,0.d0),tmp1,n_pbtot_red,
c     &                 v2_q(1,1,iq2),n_pbtot_red,(1.d0,0.d0),
c     &                 t_nu(1,1,i_nu1),np)
          enddo  !! over i_nu1
c ---------------------------------------------------------------
	      if(in==1) t_nu=conjg(t_nu)
          allocate(q_tau(np,np3,0:n_tau))
          call zgemm('n','n',np*np3,n_tau+1,nomnu,(1.d0,0.d0),t_nu,
     &               np*np3,tau_nu_from_nu_nu_1(1,0,ind_nu),
     &               n_nu2_max,(0.d0,0.d0),q_tau,np*np3)
          deallocate(t_nu)
          if(in==1) q_tau=conjg(q_tau)
c  ---------- Add static asymptotics -----------------------
          do i0_tau=0,n_tau/2
            do it=1,2
              i_tau=i0_tau
              if(it==2) i_tau=n_tau-i0_tau
              qt_pw(1:np,1:np3,iq1,in,it,i0_tau)=q_tau(:,:,i_tau)
c	          call zgemm('n','n',np,np2,np,(1.d0,0.d0),v2_q(1,1,iq1),
c     &	                 n_pbtot_red,q_pw(1,1,it,i0_tau,in),
c     &                   n_pbtot_red,(0.d0,0.d0),tmp1,n_pbtot_red)
c	          call zgemm('n','n',np,np2,np2,(1.d0,0.d0),tmp1,
c     &	                 n_pbtot_red,v2_q(1,1,iq2),n_pbtot_red,
c     &                   (1.d0,0.d0),qt_pw(1,1,iq1,in,it,i0_tau),
c     &                   n_pbtot_red)
            enddo
          enddo
          deallocate(q_tau)
        enddo  !! over in
      enddo  !! over iq1
      deallocate(tmp,tmp1,q_pw,w_nu_nu,wexa,wspl_geom,wspl_asy,wasy)
      end
