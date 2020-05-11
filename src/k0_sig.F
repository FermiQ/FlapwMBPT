      subroutine k0_sig(ind_omega,k1,k_pw,g_om_om,ps3,s2,nom,ispin,wnu)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,k1,s2,nom,ispin
      real*8, intent(in) :: wnu(n_pbtot_red,n_pbtot_red,0:n_nu,npnt_c)
      complex*16, intent(in) :: g_om_om(nbndf_bnd,nbndf_bnd,2,nom,
     &                                  npnt_c),
     &                          ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,
     &                              nqdiv_c)
	  complex*16, intent(out) :: k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c,
     &                                0:n_tau,2)
      integer :: i_omega,k,kq,nn,it,n,k0,nk,kq0,i1_omega,np,i_tau,
     &           nkq,kf0,kaq0
      real*8 :: v(3),w1,pib
      real*8, allocatable :: wexa(:,:,:),wspl_geom(:,:,:,:),
     &                       wspl_asy(:,:,:,:),wasy(:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:),tmp2(:,:,:),
     &                           w_om_om(:,:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
      n=nbndf_bnd
      np=n_pbtot_red
      allocate(tmp(n,np))
      allocate(w_om_om(n_pbtot_red,n_pbtot_red))
      allocate(wexa(n_pbtot_red,n_pbtot_red,0:n_nu_exa))
      allocate(wspl_geom(n_pbtot_red,n_pbtot_red,0:n_nu_geom+2,4))
      allocate(wspl_asy(n_pbtot_red,n_pbtot_red,n_nu_asy+1,4))
      allocate(wasy(n_pbtot_red,n_pbtot_red))
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
        call w_spl_q_red(kq,wnu,wexa,wspl_geom,wspl_asy,wasy)
        allocate(tmp1(nk,nkq,n_omega1_max)) 
        allocate(tmp2(nk,nkq,0:n_tau))   
	    do it=1,2
	      tmp1=(0.d0,0.d0)
	      do i_omega=1,nom
            if(it==1) w1=w_om_adapt_omega(i_omega,i1_omega)+pib
            if(it==2) w1=w_omega(i1_omega)
     &                  -w_om_adapt_omega(i_omega,i1_omega)
            call interp_w_nu_sq_spl(wexa,wspl_geom,wspl_asy,wasy,
     &                              w_om_om,w1,n_pbtot_red)
	        call zgemm('n','c',nk,nkq,nk,(1.d0,0.d0),
     &                 g_om_om(1,1,it,i_omega,k0),n,ps3(1,1,s2,k),np,
     &                 (0.d0,0.d0),tmp,n)
	        call zgemm('n','t',nk,nkq,nkq,(1.d0,0.d0),tmp,n,
     &                 w_om_om,n_pbtot_red,(0.d0,0.d0),
     &                 tmp1(1,1,i_omega),nk)
            if(it==1) tmp1(:,:,i_omega)=conjg(tmp1(:,:,i_omega))
          enddo  !! over i_omega
          call zgemm('n','n',nn,n_tau+1,nom,(1.d0,0.d0),tmp1,nn,
     &	             tau_om_from_omega_om(1,0,ind_omega),n_omega1_max,
     &               (0.d0,0.d0),tmp2,nn)
          if(it==1) tmp2=conjg(tmp2)
          do i_tau=0,n_tau
            k_pw(1:nk,1:nkq,k,i_tau,it)=tmp2(:,:,i_tau)
          enddo
        enddo
        deallocate(tmp1,tmp2)
      enddo   !! over k
      deallocate(tmp,w_om_om,wexa,wspl_geom,wspl_asy,wasy)
      end
