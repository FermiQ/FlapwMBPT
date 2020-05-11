      subroutine k_sig(ind_omega,k1,k_pw,g_om_om,ps3,nom,ispin,wnu,
     &                 numtn,t_pw,t_stat,v2_q,itr,sig_om)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,k1,nom,ispin,itr,numtn
      real*8, intent(in) :: wnu(n_pbtot_red,n_pbtot_red,0:n_nu,npnt_c)
      complex*16, intent(in) :: g_om_om(nbndf_bnd,nbndf_bnd,2,nom,
     &                                  npnt_c),
     &                          ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,
     &                              nqdiv_c),
     &                          t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c,
     &                               numtn,2),
     &                          t_stat(n_pbtot_red,nbndf_bnd,nqdiv_c,2,
     &                                 0:n_tau/2),
     &                          v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c)
	  complex*16, intent(inout) :: sig_om(nbndf_bnd)
	  complex*16, intent(out) :: k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c,
     &                                0:n_tau,2)
      integer :: i_omega,k,kq,nn,it,n,k0,nk,kq0,i1_omega,np,i_tau,
     &           kf0,kaq0,nk0,npb,j,i,ii
      real*8 :: v(3),w1,pib,fc
      real*8, allocatable :: wexa(:,:,:),wspl_geom(:,:,:,:),
     &                       wspl_asy(:,:,:,:),wasy(:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:),tmp2(:,:,:),
     &                           w_om_om(:,:),x(:,:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
      n=nbndf_bnd
      np=n_pbtot_red
      allocate(w_om_om(n_pbtot_red,n_pbtot_red))
      allocate(wexa(n_pbtot_red,n_pbtot_red,0:n_nu_exa))
      allocate(wspl_geom(n_pbtot_red,n_pbtot_red,0:n_nu_geom+2,4))
      allocate(wspl_asy(n_pbtot_red,n_pbtot_red,n_nu_asy+1,4))
      allocate(wasy(n_pbtot_red,n_pbtot_red))
c ------------------------------------------------------------------
      kf0=k_a_from_c(k1)
	  nk0=n_low_bnd(kf0,ispin)
      fc=(-1.d0)**itr/nqdiv_c
      k_pw=(0.d0,0.d0)
      do k=1,nqdiv_c
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
	    v=pnt_c(:,k1)-pnt_c(:,k)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K'-K
	    kq0=i_kref_c(kq)
        kaq0=k_a_from_c(kq0)
	    npb=n_pbmt_red+nplwgw_red(kaq0)
        nn=nk*npb
        call w_spl_q_red(kq,wnu,wexa,wspl_geom,wspl_asy,wasy)
        allocate(tmp1(npb,nk,n_omega1_max)) 
        allocate(tmp2(npb,nk,0:n_tau))
        allocate(x(npb,nk))
        allocate(tmp(npb,nk))
	    do it=1,2
	      tmp1=(0.d0,0.d0)
	      do i_omega=1,nom
            if(it==1) w1=w_om_adapt_omega(i_omega,i1_omega)+pib
            if(it==2) w1=w_omega(i1_omega)
     &                  -w_om_adapt_omega(i_omega,i1_omega)
            call interp_w_nu_sq_spl(wexa,wspl_geom,wspl_asy,wasy,
     &                              w_om_om,w1,n_pbtot_red)
            w_om_om=w_om_om+v2_q(:,:,kq)
            x=t_pw(1:npb,1:nk,kq,i_omega,it)
            if(itr==2) then
              call tsig_stat_tau_to_omega(t_stat,kq,npb,nk,x,i_omega,it,
     &                    omega_om_from_tau_om(0,1,1,it,ind_omega))
            endif
            call zgemm('n','n',npb,nk,npb,(1.d0,0.d0),
     &                 w_om_om,n_pbtot_red,x,npb,(0.d0,0.d0),tmp,npb)
	        call zgemm('n','t',npb,nk,nk,(1.d0,0.d0),tmp,npb,
     &	               g_om_om(1,1,it,i_omega,k0),nbndf_bnd,(0.d0,0.d0),
     &                 tmp1(1,1,i_omega),npb)
            if(it==1) tmp1(:,:,i_omega)=conjg(tmp1(:,:,i_omega))
          enddo  !! over i_omega
          call zgemm('n','n',nn,n_tau+1,nom,(1.d0,0.d0),tmp1,nn,
     &	             tau_om_from_omega_om(1,0,ind_omega),n_omega1_max,
     &               (0.d0,0.d0),tmp2,nn)
          if(it==1) tmp2=conjg(tmp2)
          do i_tau=0,n_tau
            do j=1,npb
              do i=1,nk
                k_pw(i,j,k,i_tau,it)=tmp2(j,i,i_tau)
              enddo
            enddo
          enddo
          do ii=1,nk0
            do j=1,npb
              do i=1,nk
                sig_om(ii)=sig_om(ii)+ps3(j,i,ii,k)*k_pw(i,j,k,0,it)*fc
              enddo
            enddo
          enddo
        enddo  !! over it
        deallocate(tmp1,tmp2,tmp,x)
      enddo   !! over k
      deallocate(w_om_om,wexa,wspl_geom,wspl_asy,wasy)
      end
      
      
      
      
      subroutine k_sig_stat(it,ind_tau,g_tau,ps3,v2_q,s2,k1,ispin,k_pw)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: it,ind_tau,k1,s2,ispin
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c)
      complex*16, intent(in) :: ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,
     &                              nqdiv_c),
     &                          v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c)
	  complex*16, intent(out) :: k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c)
      integer :: k,kq,nn,n,k0,nk,kq0,np,nkq,kf0,kaq0
      real*8 :: v(3)
      complex*16, allocatable :: tmp(:,:),gg(:,:)
      n=nbndf_bnd
      np=n_pbtot_red
      allocate(tmp(n,np))
c ------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
      do k=1,nqdiv_c
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
	    v=pnt_c(:,k1)-pnt_c(:,k)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K'-K
	    kq0=i_kref_c(kq)
        kaq0=k_a_from_c(kq0)
	    nkq=n_pbmt_red+nplwgw_red(kaq0)
        nn=nk*nkq
        allocate(gg(nk,nk))
        call ferm_unpack_tau(gg,g_tau(1,1,1,ind_tau,k0),nk,nk,nbndf_bnd,
     &                       it)
        call zgemm('n','c',nk,nkq,nk,(1.d0,0.d0),gg,nk,ps3(1,1,s2,k),np,
     &             (0.d0,0.d0),tmp,n)
        call zgemm('n','t',nk,nkq,nkq,(1.d0,0.d0),tmp,n,v2_q(1,1,kq),np,
     &             (0.d0,0.d0),k_pw(1,1,k),nbndf_bnd)
        deallocate(gg)
      enddo   !! over k
      deallocate(tmp)
      end
