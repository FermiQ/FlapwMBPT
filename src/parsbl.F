      subroutine parsbl(t,nt,i)
	  integer :: nt,i
      character*1 t(nt)
  99  if(t(i).ne.'_') return
      i=i+1
      if(i.gt.nt) return
      goto 99
      end
      
      
      
      
      subroutine parspi(t,nt,i,integr,ndig)
	  integer :: nt,i,integr,ndig,is
      character*1 t(nt),csym(10)
      data csym/'0','1','2','3','4','5','6','7','8','9'/
      integr=0
      ndig=0
      i=i-1
  99  i=i+1
      do is=1,10
        if(csym(is).eq.t(i)) then
          integr=10*integr+is-1
          ndig=ndig+1
          goto 99
        endif
      enddo
      end
      
      subroutine parssr(t,nt,i,r)
	  integer :: nt,i,isig,intg1,ndig,intg2
      character*1 :: t(nt)
      real*8 :: r
      isig=1
      if(t(i).eq.'-') isig=-1
      if(t(i).eq.'-'.or.t(i).eq.'+') i=i+1
      call parspi(t,nt,i,intg1,ndig)
      if(t(i).eq.'.') then
        i=i+1
        call parspi(t,nt,i,intg2,ndig)
        r=isig*(intg1+intg2*(0.1d0**ndig))
      else
        r=intg1*isig
      endif
      end
      
      subroutine parsvc(t,nt,i,v,rbas)
	  use manager_mod
	  use parallel_mod
      character*80, intent(in) :: t
      integer, intent(in) :: nt
	  integer :: i,ii,integr,ndig,m,jj
      real*8 :: v(3),rbas(3,3),fac,r3(3)
      v(1)=0.d0
      v(2)=0.d0
      v(3)=0.d0
      if(t(i:i).eq.'X') v(1)=1.d0
      if(t(i:i).eq.'Y') v(2)=1.d0
      if(t(i:i).eq.'Z') v(3)=1.d0
      if(t(i:i).eq.'X'.or.t(i:i).eq.'Y'.or.t(i:i).eq.'Z') goto 90
      if(t(i:i).eq.'D') then
        v(1)=1.d0
        v(2)=1.d0
        v(3)=1.d0
        goto 90
      endif
      if(t(i:i).eq.'S') then
	    v(1)=0.d0
	    v(2)=0.d0
	    v(3)=1.d0
        i=i+1
        call parspi(t,nt,i,integr,ndig)
        if(integr.ge.10) then
          fac=dfloat(mod(integr,10))/dfloat(integr/10)
        else
          fac=1.d0/dfloat(integr)
        endif
        do ii=1,3
          v(ii)=v(ii)*fac
        enddo
        i=i-1
        goto 90
      endif
      if(t(i:i).ne.'(') then
        if(maswrk) write(iun,*)
     &       '*** vec must be "X","Y","Z","D" or start with "("'
        call ending
      endif
      i1=i+1
      i2=i1+index(t(i1:nt),')')-2
      call decode_vector(t,i1,i2,r3)
      do m=1,3
 	    do jj=1,3
	      v(jj)=v(jj)+rbas(jj,m)*r3(m)
	    enddo
      enddo
      i=i2+1
  90  i=i+1
      end
