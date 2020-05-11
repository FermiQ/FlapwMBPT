      subroutine zerov(vatm,vmth,vinh)
***********************************************
*   set up muffin-tin zero for the potential  *
***********************************************
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      real*8, intent(inout) :: vatm(maxnsym,nsort),vmth(maxmt/nspin)
      complex*16, intent(inout) :: vinh(nplwro)
	  integer :: isort,mt,irad
	  real*8 :: sqpi4,vmt0
      sqpi4=2.d0*sqrt(pi)
c  define mt-zero
      vmt0=0.d0
      do isort=1,nsort
        mt=indmt(1,isort,1)
        vmt0=vmt0+(vmth(mt+nrad(isort))/sqpi4
     &                        -2.d0*z(isort)/smt(isort))*nhsort(isort)
      enddo
      vmt0=vmt0/dfloat(natom)
c  shift potential according to new zero
      do isort=1,nsort
        mt=indmt(1,isort,1)
        do irad=mt,mt+nrad(isort)
          vmth(irad)=vmth(irad)-vmt0*sqpi4
        enddo
        vatm(1,isort)=vatm(1,isort)-vmt0*sqpi4
      enddo
      vinh(1)=vinh(1)-vmt0
      if(maswrk) write(iun,*)' MT-zero is ',vmt0
      end
