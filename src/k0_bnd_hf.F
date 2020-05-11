      subroutine k0_bnd_hf(i_nu,iq,k0_pw,ps3,s2,g_x,ex)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_nu,iq,s2
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c,nspin),
     &                      ex(nbndf_bnd,npnt_c,nspin)
      complex*16, intent(in) :: ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,
     &                              nqdiv_c,nspin)
      complex*16, intent(out) :: k0_pw(nbndf_bnd,nbndf_bnd,nqdiv_c,
     &                                 nspin)
      integer :: k,kq,ispin,k0,nk,kq0,nkq,i,j,kf0,kfq0
      real*8 :: v(3),de
      complex*16 :: kx
c ------------------------------------------------------------------
      k0_pw=(0.d0,0.d0)
      do ispin=1,nspin
        do k=1,nqdiv_c
          k0=i_kref_c(k)
          kf0=k_a_from_c(k0)
          nk=n_low_bnd(kf0,ispin)
          v=pnt_c(:,k)-pnt_c(:,iq)
          call zone1_number(v,rb0_c,ndiv_c,kq)
          kq=index_k1_c(kq)     !! for K-Q
          kq0=i_kref_c(kq)
          kfq0=k_a_from_c(kq0)
          nkq=n_low_bnd(kfq0,ispin)
          do j=1,nkq
            do i=1,nk
              de=ex(j,kq0,ispin)-ex(i,k0,ispin)
              if(i_nu==0.and.abs(de)<1.d-8) then
                kx=-betta_t*g_x(i,n_tau,k0,ispin)
     &                     *(1.d0+g_x(i,n_tau,k0,ispin))
              else
                kx=(g_x(j,n_tau,kq0,ispin)
     &             -g_x(i,n_tau,k0,ispin))/cmplx(de,w_nu(i_nu))
              endif
              k0_pw(i,j,k,ispin)=kx*ps3(i,j,s2,k,ispin)
            enddo
          enddo
        enddo   !! over k
      enddo   !! over ispin
      end
