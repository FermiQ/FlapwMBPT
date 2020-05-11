      subroutine g_red_k_asy(g1_as,z_red)
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd,npnt,nspin)
      complex*16, intent(out) :: g1_as(nfun_red,nfun_red,nqdiv,nspin)
      integer :: ispin,ind_k,k,k0,nb,ig,k1
      complex*16, allocatable :: gk(:,:)
	allocate(gk(nfun_red,nfun_red))
	g1_as=(0.d0,0.d0)
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me_k+1)
	    k=n3_mpi_k(me_k+1)+ind_k
          nb=n_low_bnd(k,ispin)
	    call zgemm('n','c',nfun_red,nfun_red,nb,(1.d0,0.d0),
     &	             z_red(1,1,k,ispin),nfun_red,z_red(1,1,k,ispin),
     &               nfun_red,(0.d0,0.d0),gk,nfun_red)
c -------------------- G1_AS on full K_mesh -------------------------
     	    do k0=1,k_star(k)
	      k1=k_list(k0,k)
	      ig=k_group(k1)
	      call sym_g_red_k(g1_as(1,1,k1,ispin),gk,k1,ig)
          enddo
        enddo  !! over ind_k
      enddo   !! over ispin
	if(nproc_k/=1) call dgop(g1_as,2*nfun_red**2*nqdiv*nspin,'  +',
     &	                       comm_k)
	deallocate(gk)
      end