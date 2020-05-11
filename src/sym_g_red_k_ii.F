	subroutine sym_g_red_k_ii(gnew,gold,k,ig)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	use vertex_mod
	implicit none
	integer, intent(in) :: k,ig
	complex*16, intent(in) :: gold(nr_full_red,nr_full_red)
	complex*16, intent(out) :: gnew(nr_full_red,nr_full_red)
	integer :: j,j0,i,i0
	real*8 :: pi2,fas,v(3)
	complex*16 :: cf
	pi2=pi+pi
	do j=1,nr_full_red
	  j0=ip_r_red(j,ig)
	  do i=1,nr_full_red
	    i0=ip_r_red(i,ig)
	    v=tshift_r_red(:,i,ig)-tshift_r_red(:,j,ig)
 	    fas=pi2*dot_product(pnt_c(:,k),v)
	    cf=dcmplx(cos(fas),sin(fas))
	    gnew(i,j)=cf*gold(i0,j0)
	  enddo   !! over i
	enddo  !! over j
	end