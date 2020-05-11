      subroutine output_m2(ch,m2,m20)
      use atom_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      character*11, intent(in) :: ch
      real*8, intent(in) :: m2(natom),m20(natom)
	integer :: isort,iatom,ig,jatom
	real*8, allocatable :: chi_sort(:),chi0_sort(:)
      allocate(chi_sort(nsort))
      allocate(chi0_sort(nsort))
	chi_sort=0.d0
	chi0_sort=0.d0
	do isort=1,nsort
	  iatom=iat_1(isort)
	  do ig=1,ngroup
	    jatom=ip(iatom,ig)
	    chi_sort(isort)=chi_sort(isort)+m2(jatom)
	    chi0_sort(isort)=chi0_sort(isort)+m20(jatom)
	  enddo
	enddo
	chi_sort=chi_sort/ngroup
	chi0_sort=chi0_sort/ngroup
c --------------------------------------------------------------------
	if(ch=='Rsp_ZZ_MT__'.and.maswrk) write(iun,*)' FDT_ZZ:'
	if(ch=='Rsp_+-_MT__'.and.maswrk) write(iun,*)' FDT_XX+YY:'
	do isort=1,nsort
	  if(maswrk) write(iun,*)txtel(isort),' : ',chi_sort(isort)
	enddo
	if(ch=='Rsp_ZZ_MT__'.and.maswrk) write(iun,*)' FDT_ZZ_0:'
	if(ch=='Rsp_+-_MT__'.and.maswrk) write(iun,*)' FDT_XX+YY_0:'
	do isort=1,nsort
	  if(maswrk) write(iun,*)txtel(isort),' : ',chi0_sort(isort)
	enddo 
      end