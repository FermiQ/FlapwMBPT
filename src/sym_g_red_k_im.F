	subroutine sym_g_red_k_im(gnew,gold,k,ig)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	use vertex_mod
	implicit none
	integer, intent(in) :: k,ig
	complex*16, intent(in) :: gold(nr_full_red,nfun_red)
	complex*16, intent(out) :: gnew(nr_full_red,nfun_red)
	integer :: iatom,isort,ind0,ind1,j,j0
	real*8 :: pi2,fas,v(3)
	complex*16 :: cf
	complex*16, allocatable :: tmp(:,:)
	pi2=pi+pi
	allocate(tmp(nr_full_red,nfun_red))
	do j=1,nr_full_red
	  call sym_left_red(ig,gold(j,:),tmp(j,:))
	enddo
	do j=1,nr_full_red
	  j0=ip_r_red(j,ig)
	  do iatom=1,natom
	    isort=is(iatom)
	    ind0=io_lem_red(iatom)
	    ind1=ind0+lfunm_red(isort)-1
	    v=tshift_r_red(:,j,ig)-tshift(:,iatom,ig)
 	    fas=pi2*dot_product(pnt_c(:,k),v)
	    cf=dcmplx(cos(fas),sin(fas))
	    gnew(j,ind0:ind1)=cf*tmp(j0,ind0:ind1)
	  enddo   !! over iatom_c
	enddo  !! over j
	deallocate(tmp)
	end