      subroutine k_bnd_hf(ispin,i_nu,k,iq,k_pw,lambda_stat,g_x,ex)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ispin,i_nu,k,iq
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c),
     &                      ex(nbndf_bnd,npnt_c)
      complex*16, intent(in) :: lambda_stat(nbndf_bnd,nbndf_bnd)
      complex*16, intent(out) :: k_pw(nbndf_bnd,nbndf_bnd)
      integer :: kq,k0,nk,kq0,nkq,j,i,kf0,kfq0
      real*8 :: v(3),de
      complex*16 :: kx
c ------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
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
          de=ex(j,kq0)-ex(i,k0)
          if(i_nu==0.and.abs(de)<1.d-8) then
            kx=-betta_t*g_x(i,n_tau,k0)*(1.d0+g_x(i,n_tau,k0))
          else
            kx=(g_x(j,n_tau,kq0)-g_x(i,n_tau,k0))/cmplx(de,w_nu(i_nu))
          endif
          k_pw(i,j)=kx*lambda_stat(i,j)
        enddo
      enddo
      end
