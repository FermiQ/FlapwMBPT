      subroutine t_reciprocal_space_static(qtpw,v2_q,iq)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c)
      complex*16, intent(inout) :: qtpw(n_pbtot_red,n_pbtot_red,nqdiv_c,
     &                                  2,0:n_tau/2)
      integer :: nn,iq1,iq2,i_nu1,iq10,np,iq20,np2,in,ka,ka2,i,i0_tau
      real*8 :: v(3)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),t_nu(:,:,:,:),
     &                           q_pw(:,:,:,:),w_n(:,:,:)
      nn=n_pbtot_red**2
c --------------- Q(tau) ---> Q(nu') -----------------------------     
      allocate(t_nu(n_pbtot_red,n_pbtot_red,0:n_nu,2))
      allocate(tmp(n_pbtot_red,n_pbtot_red))
      allocate(tmp1(n_pbtot_red,n_pbtot_red))
      allocate(q_pw(n_pbtot_red,n_pbtot_red,2,0:n_tau/2))
      allocate(w_n(n_pbtot_red,n_pbtot_red,2))
      do iq1=1,nqdiv_c
        do i0_tau=0,n_tau/2
          do i=1,2
	      q_pw(:,:,i,i0_tau)=qtpw(:,:,iq1,i,i0_tau)
          enddo
        enddo
	    call tau_to_nu_full_kk4(t_nu,q_pw,nn)
c -----------------------------------------------------------
	    iq10=i_kref_c(iq1)
        ka=k_a_from_c(iq10)
	    np=n_pbmt_red+nplwgw_red(ka)
        v=pnt_c(:,iq1)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,iq2)  ! iq2 - index in BZ_1
        iq2=index_k1_c(iq2)  ! Now iq2 - index in BZ_0
	    iq20=i_kref_c(iq2)
        ka2=k_a_from_c(iq20)
	    np2=n_pbmt_red+nplwgw_red(ka2)
c ---------- Get T_q_nu and subtract static asymptotics ----------
        do i_nu1=0,n_nu
          call unpack_hermit(w_n(1,1,1),wnu_vrt(1,1,i_nu1,iq10),np,
     &                       n_pbtot_red,n_pbtot_red,(0.d0,0.d0),
     &                       (1.d0,0.d0))
	      call sym_w_red(iq1,iq10,w_n(1,1,1))
          w_n(:,:,1)=w_n(:,:,1)+v2_q(:,:,iq1)
          call unpack_hermit(w_n(1,1,2),wnu_vrt(1,1,i_nu1,iq20),np,
     &                       n_pbtot_red,n_pbtot_red,(0.d0,0.d0),
     &                       (1.d0,0.d0))
	      call sym_w_red(iq2,iq20,w_n(1,1,2))
          w_n(:,:,2)=w_n(:,:,2)+v2_q(:,:,iq2)
	      do in=1,2   !  -+ nu'
	        call zgemm('n','n',np,np2,np,(1.d0,0.d0),w_n(1,1,1),
     &	               n_pbtot_red,t_nu(1,1,i_nu1,in),n_pbtot_red,
     &                 (0.d0,0.d0),tmp,n_pbtot_red)
	        call zgemm('n','n',np,np2,np2,(1.d0,0.d0),tmp,n_pbtot_red,
     &	               w_n(1,1,2),n_pbtot_red,(0.d0,0.d0),
     &                 t_nu(1,1,i_nu1,in),n_pbtot_red)
          enddo
        enddo  !! over ind_nu1
c --------------- T(nu') ---> T(tau) -----------------------------
	    call nu_to_tau_kk4(t_nu,q_pw,nn)
        do i0_tau=0,n_tau/2
          do i=1,2
	        qtpw(:,:,iq1,i,i0_tau)=q_pw(:,:,i,i0_tau)
          enddo
        enddo
      enddo  !! over iq1
      deallocate(tmp,tmp1,t_nu,q_pw,w_n)
      end
