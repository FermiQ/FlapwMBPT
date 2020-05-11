      subroutine lat_sums(q,idt,lmax,dl,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
c     reduced structure constants by ewald method.
c     functions here are real hankel functions without extra factors
c     E=0  q.ne.0
      implicit none
      integer, intent(in) :: idt,lmax,khbr
      real*8, intent(in)  :: q(3)
      complex*16, intent(out) :: dl((lmax+1)**2)
      real*8 :: pi2
c ------------------------------------------------------------------
      pi2=pi+pi
c ------------------------------------------------------------------
      if(khbr==0) then
        call lat_q_sum(q,idt,lmax,dl)
        call lat_r_sum(q,dtau(1,idt),lmax,dl)
      else
        call lat_q_sum_b(q,idt,lmax,dl)
        call lat_r_sum_b(q,dtau(1,idt),lmax,dl)
      endif
      end
c --------------------------------------------
      subroutine lat_q_sum(q,idt,lmax,dl)
c     k-space part of reduced structure constants (Ewald).
c     tau in alat-units. q and rlat in tpiba-units.
c ------- We calculate only regular contribution here --------------
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: idt,lmax
      complex*16, intent(out) :: dl((lmax+1)**2)
      integer :: nlm,ilm,ir,l,m
      real*8 :: q(3),dta(3),r(3),tpi,gamma,fpibv,tpiba,qq,scalp,r2,r1,
     &          den0,sqpi4,q2a
      complex*16 :: eiphi,cs
      real*8, allocatable :: yl(:)
      dta=dtau(:,idt)
      tpi=2.d0*pi
      gamma=0.25d0/(ewald*ewald)   !!  = 1/(4*eta^2)
      fpibv=2.d0*tpi/amega   !!  = 4pi/omega
      tpiba=tpi/par      !!  = 2pi/a
      sqpi4=sqrt(4.d0*pi)
      nlm=(lmax+1)**2
      qq=dot_product(q,q)
      dl=(0.d0,0.d0)
      q2a=q2aver
      allocate(yl(nlm))
      do ir=1,nkr
        r=tpiba*(q+rlat(:,ir))
        scalp=par*dot_product(r,dta)
        eiphi=dcmplx(dcos(scalp),dsin(scalp))
        r2=dot_product(r,r)
        if(sqrt(r2)>1.d-8) then
          den0=dexp(-gamma*r2)/r2
          r1=dsqrt(r2)
          call sphharm(r(1)/r1,r(2)/r1,r(3)/r1,lmax,yl)
          ilm=0
          do l=0,lmax
            cs=eiphi*den0*r1**l
c The above line is = !k+G!^l*e^(i(k+G)delt)*e^(-(k+G)^2/4eta^2)/!k+G!^2
            do m=-l,l
              ilm=ilm+1
              dl(ilm)=dl(ilm)+yl(ilm)*cs
            enddo
          enddo
        endif
      enddo
      eiphi=(0.d0,-1.d0)*fpibv
      ilm=0
      do l=0,lmax
        eiphi=eiphi*(0.d0,1.d0)
        do m=-l,l
          ilm=ilm+1
          dl(ilm)=eiphi*dl(ilm)
        enddo
      enddo
      deallocate(yl)
c correction from Taylor expansion of the exponent ------
      if(dot_product(q,q)<1.d-8) then
        dl(1)=dl(1)-sqrt(pi)/2/amega/ewald**2
      endif
      end





      subroutine lat_r_sum(q,dta,lmax,dl)
c     adds real space part of reduced structure constants (ewald).
      use solid_mod
      use units_mod
      implicit none
      integer :: lmax,ir1,ir,l,m,ilm
      real*8 :: r(3),chi(0:40),q(3),dta(3),tpi,srpi,ta,a2,
     &        ta2,cc,r1,r2,gl,qdotr,g0,erfc_loc
      complex*16 :: cfac,dl((lmax+1)**2),zz
      real*8, allocatable :: yl(:)
      tpi=2.d0*pi
      srpi=dsqrt(pi)
      ta=2.d0*ewald     !! = 2*eta
      a2=ewald*ewald        !! = eta^2
      ta2=2.d0*a2   !! = 2*eta^2
      cc=4.d0*a2*ewald/srpi  !! = 4*eta^3/sqrt(pi)
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
        chi(0)=sqrt(pi)*erfc_loc(ewald*r1)/r1
        gl=2.d0*ewald*dexp(-a2*r2)
        do l=1,lmax
          chi(l)=((2*l-1)*chi(l-1)+gl)/r2/2.d0
          gl=a2*gl
        enddo
        qdotr=tpi*(q(1)*dlat(1,ir)+q(2)*dlat(2,ir)+q(3)*dlat(3,ir))
        cfac=dcmplx(dcos(qdotr),-dsin(qdotr))
        ilm=0
        do l=0,lmax
          zz=chi(l)*cfac*(-2*r1)**l/sqrt(pi)
          do m=-l,l
            ilm=ilm+1
            dl(ilm)=dl(ilm)+yl(ilm)*zz
          enddo
        enddo
      enddo
      deallocate(yl)
      if(ir1.eq.1) return
c --- add dl3 for diagonal structure constants ------
      g0=-ewald/pi
      dl(1)=dl(1)+g0
      end
