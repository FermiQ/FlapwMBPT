      subroutine rr_mesh_coarse
	use manager_mod
	use parallel_mod
      use sitegen_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer :: i1,i2,i3,ir
c -----  Coarse R-vectors inside parallelepiped ------------------
      allocate(rr_coarse(3,nqdiv_c))
      ir=0
      do i3=0,ndiv_c(3)-1
        do i2=0,ndiv_c(2)-1
          do i1=0,ndiv_c(1)-1
            ir=ir+1
            rr_coarse(:,ir)=i1*rbas(:,1)+i2*rbas(:,2)+i3*rbas(:,3)
          enddo
        enddo
      enddo
      end
