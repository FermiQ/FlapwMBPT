	subroutine sym_g_red_k(gnew,gold,k,ig)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	use vertex_mod
	implicit none
	integer, intent(in) :: k,ig
	complex*16, intent(in) :: gold(nfun_red,nfun_red)
	complex*16, intent(out) :: gnew(nfun_red,nfun_red)
	integer :: iatom,jatom,isort,ind0,jnd0,ind1,j,jsort,jnd1
	real*8 :: pi2,fas,v(3)
	complex*16 :: cf
	complex*16, allocatable :: tmp(:,:)
	pi2=pi+pi
	allocate(tmp(nfun_red,nfun_red))
	do j=1,nfun_red
	  call sym_left_red(ig,gold(1,j),tmp(1,j))
	enddo
	do j=1,nfun_red
	  call sym_left_red(ig,tmp(j,:),gnew(j,:))
	enddo
	do jatom=1,natom
	  jsort=is(jatom)
	  jnd0=io_lem_red(jatom)
	  jnd1=jnd0+lfunm_red(jsort)-1
	  do iatom=1,natom
	    isort=is(iatom)
	    ind0=io_lem_red(iatom)
	    ind1=ind0+lfunm_red(isort)-1
	    v=tshift(:,iatom,ig)-tshift(:,jatom,ig)
 	    fas=pi2*dot_product(pnt_c(:,k),v)
	    cf=dcmplx(cos(fas),sin(fas))
	    gnew(ind0:ind1,jnd0:jnd1)=cf*gnew(ind0:ind1,jnd0:jnd1)
	  enddo   !! over iatom_c
	enddo  !! over jatom_c
	deallocate(tmp)
	end