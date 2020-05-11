      subroutine k_real_space_mm_im(iq,k,k_pw,kcom,z_red,isp,jsp)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,isp,jsp,k
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(out) :: kcom(nbndf_bnd,nfun_red)
      integer :: k0,nk,kq,kq0,nkq,kf0,kfq0
      real*8 :: v(3)
      k0=i_kref_c(k)
      kf0=k_a_from_c(k0)
	  nk=n_low_bnd(kf0,isp)
	  v=pnt_c(:,k)-pnt_c(:,iq)
	  call zone1_number(v,rb0_c,ndiv_c,kq)
	  kq=index_k1_c(kq)     !! for K-Q
	  kq0=i_kref_c(kq)
      kfq0=k_a_from_c(kq0)
	  nkq=n_low_bnd(kfq0,jsp)
	  call zgemm('n','c',nk,nfun_red,nkq,(1.d0,0.d0),k_pw,nbndf_bnd,
     &           z_red(1,1,kq,jsp),nfun_red,(0.d0,0.d0),kcom,nbndf_bnd)
      end
