      subroutine v0_b0_finder
      use units_mod
      implicit none
      integer :: i,n,i1,i2,nn,nat
      real*8 :: v1,v2,d1,d2,v0,e0,e2,hh,v3,d3,d
      real*8, allocatable :: vv0(:),e(:)
      nn=100
      open(3,file='ini_v0b0')
      read(3,*)n,v0,i1,nat
      allocate(vv0(n))
      allocate(e(n))
      do i=1,n
        read(3,*)vv0(i),e(i)
      enddo
      close(3)
      i2=i1+1
      v1=vv0(i1)
      v2=vv0(i2)
      hh=(v2-v1)/nn
      call mnk_v0b0(n,v1,vv0,e,e0,e2,d1)
      do i=1,nn
        v2=v1+hh*i
        call mnk_v0b0(n,v2,vv0,e,e0,e2,d2)
        if(d2>=d1) then
          v3=v2-hh
          goto 1
        endif
        d1=d2
      enddo
      v3=v2
1     call mnk_v0b0(n,v3,vv0,e,e0,e2,d3)
      open(3,file='out_v0b0')
      write(3,'(a10,f12.5,1x,f12.5,a5,f12.5,a4)')' vv0 V0 = ',v3,v3*v0,
     &      ' a.u.',v3*v0*0.529177d0**3,' A^3'
      d=(0.75d0*v3*v0/pi/nat)**(1.d0/3.d0)
      write(3,*)' SWZ = ',d,' a.u.',d*0.529177d0,' A'
      write(3,*)' B0 = ',2.d0*e2/v0/6.7977d-6,' kbar'
      do i=1,n
        v1=vv0(i)-v3
        d=e0+e2*v1**2
        write(3,'(f5.3,2(1x,f15.8))')vv0(i),e(i),d
      enddo
      close(3)
      deallocate(vv0,e)
      stop
      end
      
      
      subroutine mnk_v0b0(n,v,vv0,e,e0,e2,d)
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: v,vv0(n),e(n)
      real*8, intent(out) :: e0,e2,d
      integer :: i
      real*8 :: x2,x4,ee,ex2,di
      real*8, allocatable :: x(:)
      allocate(x(n))
      x=vv0-v
      x2=0.d0
      x4=0.d0
      ee=0.d0
      ex2=0.d0
      do i=1,n
        x2=x2+x(i)**2
        x4=x4+x(i)**4
        ee=ee+e(i)
        ex2=ex2+e(i)*x(i)**2
      enddo
      di=n*x4-x2*x2
      e0=(ee*x4-ex2*x2)/di
      e2=(n*ex2-x2*ee)/di
      d=0.d0
      do i=1,n
        d=d+(e(i)-e0-e2*x(i)**2)**2
      enddo
      d=sqrt(d/n)
      deallocate(x)
      end