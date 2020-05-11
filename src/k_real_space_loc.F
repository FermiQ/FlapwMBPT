      subroutine k_real_space_loc(kcom,k_rs,z_red,isp)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: isp
      complex*16, intent(in) :: kcom(nbndf_bnd,nfun_red,nqdiv_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(out) :: k_rs(maxel_red,maxel_red,natom)
      integer :: k0,nk,kf0,iatom,isort,i0,n,k
      k_rs=(0.d0,0.d0)
      do k=1,nqdiv_c
        k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	  nk=n_low_bnd(kf0,isp)
        do iatom=1,natom
          i0=io_lem_red(iatom)
          isort=is(iatom)
          n=lfunm_red(isort)
          call zgemm('n','n',n,n,nk,(1.d0,0.d0),z_red(i0,1,k,isp),
     &               nfun_red,kcom(1,i0,k),nbndf_bnd,(1.d0,0.d0),
     &               k_rs(1,1,iatom),maxel_red)
        enddo
      enddo
      k_rs=k_rs/nqdiv_c
      end