      subroutine lat_sum_b(q,dta,lmax,dl)
c     adds real space part of reduced structure constants (ewald).
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmax
      real*8, intent(in) :: q(3),dta(3)
      complex*16, intent(out) :: dl((lmax+1)**2)
      integer :: ir1,ir,l,m,ilm
      real*8 :: r(3),tpi,r1,r2,qdotr
      complex*16 :: cfac,zz
      real*8, allocatable :: yl(:),cl(:)
      tpi=2.d0*pi
      ir1=2
      if(dta(1)**2+dta(2)**2+dta(3)**2.gt.1.d-7) ir1=1
      allocate(yl((lmax+1)**2))
      allocate(cl(0:lmax))
      dl=(0.d0,0.d0)
      do ir=ir1,nkd
        r(1)=par*(dlat(1,ir)+dta(1))
        r(2)=par*(dlat(2,ir)+dta(2))
        r(3)=par*(dlat(3,ir)+dta(3))
        r2=dot_product(r,r)
        r1=dsqrt(r2)
        call sphharm(r(1)/r1,r(2)/r1,r(3)/r1,lmax,yl)
        call c_bessel(omega_hse,r1,lmax,cl)
        qdotr=tpi*(q(1)*dlat(1,ir)+q(2)*dlat(2,ir)+q(3)*dlat(3,ir))
        cfac=dcmplx(dcos(qdotr),-dsin(qdotr))
        ilm=0
        do l=0,lmax
          zz=cl(l)*cfac
          do m=-l,l
            ilm=ilm+1
            dl(ilm)=dl(ilm)+yl(ilm)*zz
          enddo
        enddo
      enddo
      deallocate(yl,cl)
      end
