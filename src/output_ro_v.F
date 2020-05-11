      subroutine output_ro_v
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: isort,i_len,nn,isym,mt,mt1,ir,i,n
	real*8 :: sqpi4
	real*8, allocatable :: tmp(:,:)
	if(.not.maswrk) return
	sqpi4=sqrt(4.d0*pi)
	allocate(tmp(0:maxnrad,5))
	i_len=len_trim(allfile)
	open(3,file=allfile(1:i_len)//'_Electronic_density')
	do isort=1,nsort
        n=nrad(isort)
	  write(3,*)' Isort ',isort
	  nn=min(5,nsym(isort))
	  do isym=1,nn
	    if(nspin==1) then
            mt=indmt(isym,isort,1)
            tmp(0:n,isym)=ro(mt:mt+n)
	    else if(nspin==2) then
            mt=indmt(isym,isort,1)
            mt1=indmt(isym,isort,2)
            tmp(0:n,isym)=ro(mt:mt+n)+ro(mt1:mt1+n)
          endif
        enddo
        do ir=0,n
          write(3,'(6(1x,e12.5))')r(ir,isort),(tmp(ir,i),i=1,nn)
        enddo
	enddo
	close(3)
	open(3,file=allfile(1:i_len)//'_Hartree_potential')
	do isort=1,nsort
        n=nrad(isort)
	  write(3,*)' Isort ',isort
	  nn=min(5,nsym(isort))
	  do isym=1,nn
          mt=indmt(isym,isort,1)
          tmp(0:n,isym)=v_mt_h(mt:mt+n)
        enddo
        do ir=0,n
          write(3,'(6(1x,e12.5))')r(ir,isort),(tmp(ir,i),i=1,nn)
        enddo
	enddo
	close(3)
	open(3,file=allfile(1:i_len)//'_Exchange_Correlation_Potential')
	do isort=1,nsort
        n=nrad(isort)
	  write(3,*)' Isort ',isort
	  nn=min(5,nsym(isort))
	  do isym=1,nn
          mt=indmt(isym,isort,1)
          tmp(0:n,isym)=v_mt_xc(mt:mt+n)
        enddo
        do ir=0,n
          write(3,'(6(1x,e12.5))')r(ir,isort),(tmp(ir,i),i=1,nn)
        enddo
	enddo
	close(3)
	deallocate(tmp)
      end