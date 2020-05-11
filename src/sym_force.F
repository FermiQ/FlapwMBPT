      subroutine sym_force(f)
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(inout) :: f(3,natom)
	integer :: ig,iatom,iatom0
      real*8 :: v(3),v1(3)
	real*8, allocatable :: f0(:,:)
      allocate(f0(3,natom))
      f0=f
      f=0.d0
      do iatom=1,natom
        do ig=1,ngroup
          iatom0=ip(iatom,ig)
          v=f0(:,iatom0)
          call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),2)
          f(:,iatom)=f(:,iatom)+v1
        enddo
      enddo
      f=f/ngroup
      deallocate(f0)
      end