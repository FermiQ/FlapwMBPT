      subroutine vertex_loc_to_bnd(vx,lambda_tau,iq,z_red,isp,jsp)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: iq,isp,jsp
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          vx(maxel_red,maxel_red,natom)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        nqdiv_c)
      integer :: k0,nk,kq,kq0,k,nkq,kf0,kfq0,iatom,isort,n,i0
      real*8 :: v(3)
      complex*16, allocatable :: tt(:,:)
	allocate(tt(maxel_red,nbndf_bnd))
c -------- Transform Vrt(tau;nu) ---- > Vrt(omega;nu) ------------
	do k=1,nqdiv_c
	  k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	  nk=n_low_bnd(kf0,isp)
	  v=pnt_c(:,k)-pnt_c(:,iq)
	  call zone1_number(v,rb0_c,ndiv_c,kq)
	  kq=index_k1_c(kq)     !! for K-Q
	  kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
	  nkq=n_low_bnd(kfq0,jsp)
        do iatom=1,natom
          i0=io_lem_red(iatom)
          isort=is(iatom)
          n=lfunm_red(isort)
          call zgemm('n','n',n,nkq,n,(1.d0,0.d0),vx(1,1,iatom),
     &               maxel_red,z_red(i0,1,kq,jsp),nfun_red,(0.d0,0.d0),
     &               tt,maxel_red)
          call zgemm('c','n',nk,nkq,n,(1.d0,0.d0),z_red(i0,1,k,isp),
     &	           nfun_red,tt,maxel_red,(1.d0,0.d0),lambda_tau(1,1,k),
     &               nbndf_bnd)
        enddo   !! over iatom
      enddo   !! over k
      deallocate(tt)
      end