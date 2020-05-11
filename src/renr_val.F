      subroutine renr_val(rhom,rhoi)
	  use atom_mod
	  use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      real*8, intent(inout) :: rhom(maxmt)
      complex*16, intent(inout) :: rhoi(nplwro,nspin)
      integer :: ispin,isort,mt,n,i
      real*8 :: zvtot,anorm,f,dqdall
      real*8, allocatable :: tm(:)
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      zvtot=0.d0
      do ispin=1,nspin
        do isort=1,nsort
          n=nrad(isort)
          f=sqrt(4.d0*pi)*nhsort(isort)
          mt=indmt(1,isort,ispin)
          allocate(tm(0:n))
          do i=0,n
            tm(i)=r(i,isort)**2*rhom(mt+i)*dr(i,isort)
          enddo
          zvtot=zvtot+f*dqdall(h(isort),tm,n)
          deallocate(tm)
        enddo   !!! over isort
        zvtot=zvtot+real(rhoi(1,ispin))*amegaint
        zvtot=zvtot+amega*dot_product(real(rhoi(2:nplwro,ispin)),
     &                                sovr(2:nplwro))
        if(complex_ro) zvtot=zvtot-amega
     &                      *dot_product(imag(rhoi(2:nplwro,ispin)),
     &	                                 sovi(2:nplwro))
      enddo   !!! over ispin
*  renormalize valence density
      anorm=nelec/zvtot
      if(abs(anorm-1.d0).gt.0.5d0) then
        if(maswrk) write(iun,*)' RENROV : Anorm = ',anorm
      endif
      do ispin=1,nspin
        do isort=1,nsort
          n=nrad(isort)
          mt=indmt(1,isort,ispin)
          do i=0,n
            rhom(mt+i)=rhom(mt+i)*anorm
          enddo
        enddo   !!! over isort
      enddo   !!! over ispin
      rhoi=anorm*rhoi
      if(maswrk) then
        write(iun,1010)nelec
        write(iun,1020)zvtot
        write(iun,1030)anorm
      endif
1010  format(1x,'valence charge in whole  elementary cell  must be ',
     &  g14.7)
1020  format(1x,'valence charge  found is ',g14.7)
1030  format(1x,'renormalization coefficient of the val.density is ',
     &  g14.7)
      end
