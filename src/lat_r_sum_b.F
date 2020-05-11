      subroutine lat_r_sum_b(q,dta,lmax,dl)
c     adds real space part of reduced structure constants (ewald).
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmax
      real*8, intent(in) :: q(3),dta(3)
      integer :: ir1,ir,l,m,ilm
      real*8 :: r(3),chi(-1:40),tpi,srpi,ta,a2,a,b,c,ta2,r1,r2,gl,
     &          qdotr,g0,erfc_loc
      complex*16 :: cfac,dl((lmax+1)**2),zz
      real*8, allocatable :: yl(:)
      tpi=2.d0*pi
      srpi=dsqrt(pi)
      ta=2.d0*ewald     !! = 2*eta
      a2=ewald*ewald        !! = eta^2
      ta2=2.d0*a2   !! = 2*eta^2
      ir1=2
      if(dta(1)**2+dta(2)**2+dta(3)**2.gt.1.d-7) ir1=1
      allocate(yl((lmax+1)**2))
      do ir=ir1,nkd
        r(1)=par*(dlat(1,ir)+dta(1))
        r(2)=par*(dlat(2,ir)+dta(2))
        r(3)=par*(dlat(3,ir)+dta(3))
        r2=dot_product(r,r)
        r1=dsqrt(r2)
        call sphharm(r(1)/r1,r(2)/r1,r(3)/r1,lmax,yl)
c ------ make the xi's from 0 to lmax ------
        a=exp(-omega_hse*r1)
        b=erfc_loc(omega_hse/ta-ewald*r1)
        c=erfc_loc(omega_hse/ta+ewald*r1)
        chi(0)=(a-0.5d0*(b*a-c/a))/r1
        chi(-1)=(a-0.5d0*(b*a+c/a))/omega_hse
        gl=ta*exp(-omega_hse**2/4.d0/a2-a2*r2)/srpi
        do l=0,lmax-1
          chi(l+1)=((2*l+1)*chi(l)+omega_hse**2*chi(l-1)+gl)/r2
          gl=ta2*gl
        enddo
        qdotr=tpi*(q(1)*dlat(1,ir)+q(2)*dlat(2,ir)+q(3)*dlat(3,ir))
        cfac=dcmplx(dcos(qdotr),-dsin(qdotr))
        ilm=0
        do l=0,lmax
          zz=chi(l)*cfac*r1**l
          do m=-l,l
            ilm=ilm+1
            dl(ilm)=dl(ilm)+yl(ilm)*zz
          enddo
        enddo
      enddo
      deallocate(yl)
      if(ir1.eq.1) return
c --- add dl3 for diagonal structure constants ------
      g0=-(ewald*exp(-omega_hse**2/4.d0/a2)
     &    -omega_hse*srpi*erfc_loc(omega_hse/ta)/2.d0)/pi
      dl(1)=dl(1)+g0
      end
