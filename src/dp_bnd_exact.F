      subroutine dp_bnd_exact(nom,i_nu,k,iq,lambda_dyn,lambda_stat,
     &                        ps3,pol,ispin,g_x,ex,trn,key,s2)
c     key = 0 for full P = P0 + dP
c     key = 1 for dP only
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nom,i_nu,k,iq,ispin,key,s2
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c),
     &                      ex(nbndf_bnd,npnt_c)
      complex*16, intent(in) :: lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2),
     &                          lambda_stat(nbndf_bnd,nbndf_bnd),
     &                          ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red),
     &                          trn(n_omega1_max)
      complex*16, intent(inout) :: pol(n_pbtot_red)
      integer :: i_omega,kq,in,k0,nk,kq0,nkq,ii,i,j,kf0,kfq0,iqa
      real*8 :: v(3),om,omnu,de
      complex*16, allocatable :: vx(:,:),tt(:,:),kx(:,:)
      iqa=k_a_from_c(iq)
c ------------------------------------------------------------------
      v=pnt_c(:,k)-pnt_c(:,iq)
      call zone1_number(v,rb0_c,ndiv_c,kq)
      kq=index_k1_c(kq)     !! for K-Q
      k0=i_kref_c(k)
      kf0=k_a_from_c(k0)
      nk=n_low_bnd(kf0,ispin)
      kq0=i_kref_c(kq)
      kfq0=k_a_from_c(kq0)
      nkq=n_low_bnd(kfq0,ispin)
      allocate(vx(nk,nkq))
      allocate(kx(nk,nkq))
c -------------------------------------------------------------------
      allocate(tt(nk,nkq))
      tt=(0.d0,0.d0)
      do in=1,2
        do i_omega=1,nom
          om=w_om_adapt_nu(i_omega,i_nu)
          omnu=om-w_nu(i_nu)
          vx=lambda_dyn(1:nk,1:nkq,i_omega,in)
          if(in==1) then
            do j=1,nkq
              do i=1,nk
                kx(i,j)=-vx(i,j)/dcmplx(-ex(i,k0),-omnu)
     &                         /dcmplx(-ex(j,kq0),-om)
              enddo
            enddo
            tt=tt+kx*conjg(trn(i_omega))
          else if(in==2) then
            do j=1,nkq
              do i=1,nk
                kx(i,j)=-vx(i,j)/dcmplx(-ex(i,k0),om)
     &                         /dcmplx(-ex(j,kq0),omnu)
              enddo
            enddo
            tt=tt+kx*trn(i_omega)
          endif
        enddo  !! over i_omega
      enddo
c --------- Add the static part ------------------------------------
      vx=lambda_stat(1:nk,1:nkq)
      if(key==0) vx=vx+ps3(1:nk,1:nkq,s2)
      do j=1,nkq
        do i=1,nk
          de=ex(j,kq0)-ex(i,k0)
          if(i_nu==0.and.abs(de)<1.d-8) then
            kx(i,j)=betta_t*vx(i,j)*g_x(i,n_tau,k0)
     &                             *(1.d0+g_x(i,n_tau,k0))
          else
            kx(i,j)=vx(i,j)*(g_x(i,n_tau,k0)-g_x(j,n_tau,kq0))
     &                 /dcmplx(de,w_nu(i_nu))
          endif
        enddo
      enddo
      tt=tt+kx
      do ii=1,n_pbmt_red+nplwgw_red(iqa)
        do j=1,nkq
          do i=1,nk
            pol(ii)=pol(ii)+conjg(ps3(i,j,ii))*tt(i,j)/nqdiv_c
          enddo
        enddo
      enddo
      deallocate(vx,tt,kx)
      end
