      subroutine vertex_wk_sing(iq,i0_tau,isp,jsp,ll,k_pw,lambda_tau)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq,i0_tau,isp,jsp,ll
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,ll,
     &                               nqdiv_c)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        nqdiv_c,ll,2)
      integer :: il,it,i_tau,k,k0,ka0,nk,kq,kq0,kqa0,nkq,j,i
      real*8 :: w,v(3)
      w=ws_head(i0_tau)*q2aver_c/amega/nqdiv_c
      do il=1,ll
        do it=1,2
          if(it==1) i_tau=i0_tau
          if(it==2) i_tau=n_tau-i0_tau
          do k=1,nqdiv_c
            k0=i_kref_c(k)
            ka0=k_a_from_c(k0)
            nk=n_low_bnd(ka0,isp)
            v=pnt_c(:,k)-pnt_c(:,iq)
            call zone1_number(v,rb0_c,ndiv_c,kq)
            kq=index_k1_c(kq)     !! for K-Q
            kq0=i_kref_c(kq)
            kqa0=k_a_from_c(kq0)
            nkq=n_low_bnd(kqa0,jsp)
            do j=1,nkq
              do i=1,nk
                lambda_tau(i,j,k,il,it)=lambda_tau(i,j,k,il,it)
     &                                 +w*k_pw(i,j,i_tau,il,k)
              enddo
            enddo
          enddo   !! over k
        enddo  !! over it
      enddo   !! over il
      end
