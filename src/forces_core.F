      subroutine forces_core(f)
	use atom_mod
	use manager_mod
      use models_mod
	use parallel_mod
	use units_mod
	use solid_mod
      use vertex_mod
      implicit none
      real*8, intent(inout) :: f(3,nsort)
      integer :: isort,m,lm,isym,mt,ir,ispin,i,mt0
      real*8 :: dqdall,fc,amro,sm,dv
      real*8, allocatable :: work(:),v(:)
      allocate(work(0:maxnrad))
      allocate(v(0:maxnrad))
      fc=1.d0/sqrt(3.d0)
      do isort=1,nsort
        if(lmpb(isort)==0) cycle
        sm=smt(isort)
        do m=-1,1
          lm=3+m
          amro=0.d0
          if(sym(lm,isort))then
            isym=lmsym(lm,isort)
            mt0=indmt(isym,isort,1)
            do ispin=1,nspin
              mt=indmt(isym,isort,ispin)
              do ir=0,nrad(isort)
                v(ir)=v_mt_h(mt0+ir)+v_mt_xc(mt+ir)
              enddo
              work(0)=0.d0
              do ir=1,nrad(isort)
                if(ir<nrad(isort)) dv=(v(ir+1)-v(ir-1))/2.d0/h(isort)
                if(ir==nrad(isort)) dv=(v(ir)-v(ir-1))/h(isort)
                dv=dv/dr(ir,isort)
                work(ir)=dr(ir,isort)*r(ir,isort)
     &                               *ro_core(ir,isort,ispin)
                work(ir)=work(ir)*(r(ir,isort)*dv+2.d0*v(ir))
              enddo
              amro=amro-fc*dqdall(h(isort),work,nrad(isort))
            enddo   !! over ispin
          endif
          if(m==-1) f(2,isort)=f(2,isort)+amro
          if(m==0) f(3,isort)=f(3,isort)+amro
          if(m==1) f(1,isort)=f(1,isort)+amro
        enddo
      enddo
      deallocate(work,v)
      if(maswrk) then
        write(iun,*)' Forces after CORE'
        do isort=1,nsort
          write(iun,'(i5,3(1x,f12.5))')isort,(f(i,isort),i=1,3)
        enddo
      endif
      end