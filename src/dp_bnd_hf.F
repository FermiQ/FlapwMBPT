      subroutine dp_bnd_hf(i_nu,k,iq,lambda_stat,ps3,pol,ispin,g_x,ex,
     &                     n0)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_nu,k,iq,ispin,n0
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c),
     &                      ex(nbndf_bnd,npnt_c)
      complex*16, intent(in) :: lambda_stat(nbndf_bnd,nbndf_bnd),
     &                          ps3(nbndf_bnd,nbndf_bnd,n0)
      complex*16, intent(inout) :: pol(n0)
      integer :: kq,k0,nk,kq0,nkq,ii,i,j,kf0,kfq0,iqa
      real*8 :: v(3),de
      complex*16, allocatable :: kx(:,:)
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
      allocate(kx(nk,nkq))
      do j=1,nkq
        do i=1,nk
          de=ex(j,kq0)-ex(i,k0)
          if(i_nu==0.and.abs(de)<1.d-8) then
            kx(i,j)=betta_t*g_x(i,n_tau,k0)*(1.d0+g_x(i,n_tau,k0))
          else
            kx(i,j)=(g_x(i,n_tau,k0)-g_x(j,n_tau,kq0))
     &             /cmplx(de,w_nu(i_nu))
          endif
          kx(i,j)=kx(i,j)*lambda_stat(i,j)
        enddo
      enddo
      do ii=1,n_pbmt_red+nplwgw_red(iqa)
        do j=1,nkq
          do i=1,nk
            pol(ii)=pol(ii)+conjg(ps3(i,j,ii))*kx(i,j)/nqdiv_c
          enddo
        enddo
      enddo
      deallocate(kx)
      end
