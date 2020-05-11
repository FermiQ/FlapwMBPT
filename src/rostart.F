      subroutine rostart(rofi,roatom)
	  use atom_mod
	  use manager_mod
	  use heg_mod
	  use solid_mod
	  use units_mod
      implicit none
      real*8, intent(in) :: roatom(0:nrmax,nsort),rofi(0:nrmax,nsort)
	  integer :: isort,mt,mt1,nsitmax,isym,irad,istar
	  real*8 :: const,radius,zvtot,zsum,dqdall,charge_int
	  real*8, allocatable :: work(:),expf(:),pref(:)
	  allocate(work(0:maxnrad))
	  if(iheg) then
	    const=ro_el_gas*sqrt(4.d0*pi)/nspin
	    ro=0.d0
	    ro_core=0.d0
	    z=0.d0
        do isort=1,nsort
          mt=indmt(1,isort,1)
          mt1=indmt(1,isort,nspin)
          ro(mt:mt+nrad(isort))=const
          ro(mt1:mt1+nrad(isort))=const
        enddo   !!! over isort
	    rointr(1,:)=ro_el_gas/nspin
	    rointr(2:,:)=0.d0
	  else
        radius=20.d0
        nsitmax=2400
	    ro=0.d0
        call mathpaw(nsitmax,radius,rofi,roatom)
        zvtot=0.d0
        zsum=0.d0
        do isort=1,nsort
          do isym=1,nsym(isort)
            mt=indmt(isym,isort,1)
            mt1=indmt(isym,isort,nspin)
            do irad=0,nrad(isort)
              ro(mt+irad)=ro(mt+irad)/nspin
              ro(mt1+irad)=ro(mt+irad)
            enddo
          enddo   !!! over isym
          mt=indmt(1,isort,1)
          const=dfloat(nspin)/sqrt(4.d0*pi)
          do irad=0,nrad(isort)
            work(irad)=ro(mt+irad)*weight(irad,isort)*const
          enddo
          zvtot=zvtot+dqdall(h(isort),work,nrad(isort))*nhsort(isort)
          zsum=zsum+z(isort)*nhsort(isort)
        enddo   !!! over isort
        call timel('******* matcor finished ************')
c     The interstitial region &&&&&&&&&&&&&&&&&&&&&&&&&&&
        charge_int=zsum-zvtot
	    rointr=(0.d0,0.d0)
        allocate(expf(nsort))
        allocate(pref(nsort))
        call PLWRHO(radius,rofi,roatom,rointr(1,1),expf,pref)
        deallocate(expf,pref)
        if(nspin.eq.2) then
          do istar=1,nplwro
            rointr(istar,1)=0.5d0*rointr(istar,1)
            rointr(istar,2)=rointr(istar,1)
          enddo
        endif
	  endif
	  deallocate(work)
	  if(magn==2) then
	    spmt=0.d0
	    spintr=(0.d0,0.d0)
	  endif
      end
