      subroutine sig_1_semi_mt(t_pw,g_tau,wnuq,ps3,s2,k1,ispin,sig_semi,
     &                         nb)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: s2,k1,ispin,nb
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c)
      complex*16, intent(in) :: wnuq(n_pbtot_red,n_pbtot_red,ndim3_nu,
     &                               nqdiv_c),
     &                          ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,
     &                              nqdiv_c),
     &                          t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c,2,
     &                               ndim3_tau)
      complex*16, intent(inout) :: sig_semi(nbndf_bnd,nbndf_bnd,2,
     &                                      ndim3_tau)
      integer :: iq,iq0,iqa0,npb,kq,kq0,kfq0,nkq,it,
     &           ind_nu,ind_tau,ii,j,i
      real*8 :: v(3)
      complex*16, allocatable :: t_nu(:,:,:,:,:),wt(:,:,:,:),
     &                           wtt(:,:,:,:),gg(:,:),tt(:,:)
      allocate(t_nu(n_pbtot_red,nbndf_bnd,nqdiv_c,ndim3_nu,2))
      call tau_to_nu_full(t_nu,t_pw,n_pbtot_red*nbndf_bnd*nqdiv_c)
      do iq=1,nqdiv_c
        iq0=i_kref_c(iq)
        iqa0=k_a_from_c(iq0)
        npb=n_pbmt_red+nplwgw_red(iqa0)
	    v=pnt_c(:,k1)-pnt_c(:,iq)
    	call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K'-Q
	    kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kfq0,ispin)
        allocate(wt(npb,nkq,ndim3_nu,2))
        do it=1,2
          do ind_nu=1,ndim3_nu
	        call zgemm('n','n',npb,nkq,npb,(1.d0,0.d0),
     &                 wnuq(1,1,ind_nu,iq),n_pbtot_red,
     &                 t_nu(1,1,iq,ind_nu,it),n_pbtot_red,
     &                 (0.d0,0.d0),wt(1,1,ind_nu,it),npb)
          enddo
        enddo
        allocate(wtt(npb,nkq,2,ndim3_tau))
        call nu_to_tau_cc(wt,wtt,npb*nkq)
        deallocate(wt)
        allocate(gg(nkq,nkq))
        allocate(tt(npb,nkq))
        do ind_tau=1,ndim3_tau
          do it=1,2
            call ferm_unpack_tau(gg,g_tau(1,1,1,ind_tau,kq0),nkq,nkq,
     &                           nbndf_bnd,it)
	        call zgemm('n','t',npb,nkq,nkq,(1.d0,0.d0),
     &                 wtt(1,1,it,ind_tau),npb,gg,nkq,(0.d0,0.d0),tt,
     &                 npb)
            do ii=1,nb
              do j=1,nkq
                do i=1,npb
                  sig_semi(ii,s2,it,ind_tau)=sig_semi(ii,s2,it,ind_tau)
     &              +ps3(i,j,ii,kq)*tt(i,j)/nqdiv_c
                enddo
              enddo
            enddo
          enddo
        enddo
        deallocate(wtt,gg,tt)
      enddo
      deallocate(t_nu)
      end
