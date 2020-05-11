      subroutine k_real_space_mm(k,kcom,k_rs,z_red,isp)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: isp,k
      complex*16, intent(in) :: kcom(nbndf_bnd,nfun_red),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(out) :: k_rs(nfun_red,nfun_red)
      integer :: k0,nk,kf0
      k0=i_kref_c(k)
      kf0=k_a_from_c(k0)
	  nk=n_low_bnd(kf0,isp)
      call zgemm('n','n',nfun_red,nfun_red,nk,(1.d0,0.d0),
     &           z_red(1,1,k,isp),nfun_red,kcom,nbndf_bnd,(0.d0,0.d0),
     &           k_rs,nfun_red)
      end
