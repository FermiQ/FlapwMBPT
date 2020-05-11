      subroutine tsig_end_dyn(nomnu,ind_omega,k1,iq,t_pw,g_nu_om,ps3,
     &                        sig_om,ispin,nk,wexa,wspl_geom,wspl_asy,
     &                        wasy,numtn)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: nomnu,ind_omega,k1,iq,ispin,nk,numtn
      real*8, intent(in) :: wexa(n_pbtot_red,n_pbtot_red,0:n_nu_exa),
     &                      wspl_geom(n_pbtot_red,n_pbtot_red,
     &                                0:n_nu_geom+2,4),
     &                      wspl_asy(n_pbtot_red,n_pbtot_red,n_nu_asy+1,
     &                               4),
     &                      wasy(n_pbtot_red,n_pbtot_red)
      complex*16, intent(in) :: g_nu_om(nbndf_bnd,nbndf_bnd,2,nomnu,
     &                                  npnt_c),
     &                          ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,
     &                              nqdiv_c),
     &                          t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c,
     &                               numtn,2)
	  complex*16, intent(inout) :: sig_om(nbndf_bnd)
      integer :: i_nu,kq,it,kq0,nkq,ii,i,j,kfq0,iq0,iqa0,npb,i1_omega
      real*8 :: v(3),pib,w1
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tt(:,:),
     &                           w_nu_om(:,:)
      pib=pi/betta_t
      i1_omega=me_t*ndim3_omega+ind_omega-1
	  v=pnt_c(:,k1)-pnt_c(:,iq)
	  call zone1_number(v,rb0_c,ndiv_c,kq)
	  kq=index_k1_c(kq)     !! for K'-Q
	  kq0=i_kref_c(kq)
      kfq0=k_a_from_c(kq0)
	  nkq=n_low_bnd(kfq0,ispin)
      iq0=i_kref_c(iq)
      iqa0=k_a_from_c(iq0)
      npb=n_pbmt_red+nplwgw_red(iqa0)
      allocate(tmp(npb,nkq))
      allocate(tmp1(npb,nkq))
      allocate(w_nu_om(n_pbtot_red,n_pbtot_red))
c -------------------------------------------------------------------
      allocate(tt(npb,nkq))
      tt=(0.d0,0.d0)     
	  do it=1,2
	    do i_nu=1,nomnu
          if(it==1) w1=w_omega(i1_omega)
     &                -w_nu_adapt_omega(i_nu,i1_omega)-pib
          if(it==2) w1=w_nu_adapt_omega(i_nu,i1_omega)
          call interp_w_nu_sq_spl(wexa,wspl_geom,wspl_asy,wasy,
     &                            w_nu_om,w1,n_pbtot_red) 
          call zgemm('n','n',npb,nkq,npb,(1.d0,0.d0),
     &               w_nu_om,n_pbtot_red,t_pw(1,1,iq,i_nu,it),
     &               n_pbtot_red,(0.d0,0.d0),tmp,npb)
	      call zgemm('n','t',npb,nkq,nkq,(1.d0,0.d0),tmp,npb,
     &	             g_nu_om(1,1,it,i_nu,kq0),nbndf_bnd,(0.d0,0.d0),
     &               tmp1,npb)
          tt=tt+tmp1*tau_om_from_nu_om_1(i_nu,0,ind_omega)
        enddo  !! over i_nu
      enddo
      deallocate(tmp,tmp1,w_nu_om)
      do ii=1,nk
        do j=1,nkq
          do i=1,npb
            sig_om(ii)=sig_om(ii)+ps3(i,j,ii,kq)*tt(i,j)/nqdiv_c
          enddo
        enddo
      enddo
      deallocate(tt)
      end
      
      
      
      subroutine tsig_end_stat(k1,iq,t_pw,g_tau,v2_q,ps3,sig,ispin,nk,
     &                         it,ind_tau)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: k1,iq,ispin,nk,it,ind_tau
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c)
      complex*16, intent(in) :: ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,
     &                              nqdiv_c),
     &                          v2_q(n_pbtot_red,n_pbtot_red),
     &                          t_pw(n_pbtot_red,nbndf_bnd)
	  complex*16, intent(inout) :: sig(nbndf_bnd)
      integer :: kq,kq0,nkq,ii,i,j,kfq0,iq0,iqa0,npb
      real*8 :: v(3)
      complex*16, allocatable :: tmp(:,:),tt(:,:),gt(:,:)
	  v=pnt_c(:,k1)-pnt_c(:,iq)
	  call zone1_number(v,rb0_c,ndiv_c,kq)
	  kq=index_k1_c(kq)     !! for K'-Q
	  kq0=i_kref_c(kq)
      kfq0=k_a_from_c(kq0)
	  nkq=n_low_bnd(kfq0,ispin)
      iq0=i_kref_c(iq)
      iqa0=k_a_from_c(iq0)
      npb=n_pbmt_red+nplwgw_red(iqa0)
      allocate(tmp(npb,nkq))
      allocate(gt(nkq,nkq))
      call ferm_unpack_tau(gt,g_tau(1,1,1,ind_tau,kq0),nkq,nkq,
     &                     nbndf_bnd,it)
c -------------------------------------------------------------------
      allocate(tt(npb,nkq))
      call zgemm('n','n',npb,nkq,npb,(1.d0,0.d0),v2_q,n_pbtot_red,t_pw,
     &           n_pbtot_red,(0.d0,0.d0),tmp,npb)
	  call zgemm('n','t',npb,nkq,nkq,(1.d0,0.d0),tmp,npb,gt,nkq,
     &           (0.d0,0.d0),tt,npb)
      deallocate(tmp,gt)
      do ii=1,nk
        do j=1,nkq
          do i=1,npb
            sig(ii)=sig(ii)+ps3(i,j,ii,kq)*tt(i,j)/nqdiv_c
          enddo
        enddo
      enddo
      deallocate(tt)
      end
