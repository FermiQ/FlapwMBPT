      subroutine forces_hf(f)
	use atom_mod
	use manager_mod
      use models_mod
	use parallel_mod
	use units_mod
	use solid_mod
      use vertex_mod
      implicit none
      real*8, intent(inout) :: f(3,nsort)
      integer :: isort,m,lm,isym,mt,mt1,irad,i
      real*8 :: dqdall,rho,fc,fc1,amro,sm
      real*8, allocatable :: work(:)
      allocate(work(0:maxnrad))
      fc=8.d0*pi/3.d0
      fc1=sqrt(3.d0/4.d0/pi)
      do isort=1,nsort
        if(lmpb(isort)==0) cycle
        sm=smt(isort)
        do m=-1,1
          lm=3+m
          amro=0.d0
          if(sym(lm,isort))then
            isym=lmsym(lm,isort)
            mt=indmt(isym,isort,1)
            mt1=indmt(isym,isort,nspin)
            do irad=0,nrad(isort)
              rho=(ro(mt+irad)+ro(mt1+irad))*nspin/2.d0
              work(irad)=rho*(1.d0-(r(irad,isort)/sm)**3)*dr(irad,isort)
            enddo
            amro=fc*dqdall(h(isort),work,nrad(isort))
     &          +vatmt(isym,isort)/sm
          endif
          if(m==-1) f(2,isort)=f(2,isort)+z(isort)*fc1*amro
          if(m==0) f(3,isort)=f(3,isort)+z(isort)*fc1*amro
          if(m==1) f(1,isort)=f(1,isort)+z(isort)*fc1*amro
        enddo
      enddo
      deallocate(work)
      if(maswrk) then
        write(iun,*)' Forces after HF'
        do isort=1,nsort
          write(iun,'(i5,3(1x,f12.5))')isort,(f(i,isort),i=1,3)
        enddo
      endif
      end