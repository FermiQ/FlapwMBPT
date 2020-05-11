      subroutine pbe_x_sr(rho,absg,t,alrho,ex,vx)
******************************************************************
*  Returns GGA96 correction to the exchange energy/potential.    *
*  for non-spin polarized system. For both spins call first      *
*  for rho=2*rho_up, then for rho=2*rho_dn.                      *
******************************************************************
      use manager_mod
      use units_mod
      IMPLICIT none
      real*8, intent(in) :: rho,absg,t,alrho
      real*8, intent(out) :: ex,vx
      real*8 :: kap,mu,kf,s,a,s2,msk,jj,jja,jjaa,a2,c,v1,v2,v3,b,hx,
     &          hxs,fx
C  GGA constants
        DATA kap/0.804000000000000000D+00/  !
     &      ,mu/0.219516451220895831D+00/  !  betta*pi^2/3
      ex=0.D0
      vx=0.D0
C  Check case rho=0
      IF(rho.LE.1.D-14) RETURN
      kf=(3.D0*PI*PI*RHO)**(1.D0/3.D0)
      s=ABSG/(2.D0*kf*RHO)
      s2=s*s
      msk=mu*s2/kap
      fx=1.d0+kap-kap/(1.d0+msk)
      a=omega_hse*sqrt(fx)/2.d0/kf
      a2=a*a
      jj=1.d0-2.d0*a2/3.d0-8.d0*a*atan(1.d0/a)/3.d0
     &       +2.d0*a2*(a2+3.d0)*log(1.d0+1.d0/a2)/3.d0
      jja=-4.d0*a/3.d0-8.d0*atan(1.d0/a)/3.d0
     &                -4.d0*(a2*a+a)/(1.d0+a2)/3.d0
     &                +4.d0*(2*a2*a+3*a)*log(1.d0+1.d0/a2)/3.d0
      jjaa=-8.d0+4.d0*(2*a+1.d0)*log(1.d0+1.d0/a2)
      c=-1.5d0*(3.d0/pi)**(1.d0/3.d0)
      ex=c*rho**(1.d0/3.d0)*fx*jj
      b=2.d0*(3.d0*pi*pi)**(1.d0/3.d0)
      hx=2.d0*mu/(1.d0+msk)**2
      hxs=-8.d0*mu*mu*s/kap/(1.d0+msk)**3
      v1=4.d0*rho**(1.d0/3.d0)*fx/3.d0-alrho*hx/b**2/rho**(4.d0/3.d0)
     &  +(4.d0*rho**(1.d0/3.d0)*s**3/3.d0-t/b**3/rho**(8.d0/3.d0))*hxs
      v2=-omega_hse*fx**1.5d0/3.d0/b
     &   +(omega_hse*s2/2.d0/b
     &    -omega_hse*alrho/2.d0/b**3/rho**(5.d0/3.d0))*sqrt(fx)*hx
     &   +(2.d0*omega_hse*s2*s/3.d0/b-omega_hse*t/2.d0/b**4/rho**3)
     &   *sqrt(fx)*hxs
     &   +(omega_hse*s2*s2/b-3.d0*omega_hse*s*t/4.d0/b**4/rho**3)*hx**2
     &    /sqrt(fx)
      v3=omega_hse**2*s2*fx*hx/6.d0/b**2/rho**(1.d0/3.d0)
     &  +omega_hse**2*(s2*s2/3.d0/b**2/rho**(1.d0/3.d0)
     &                -s*t/4.d0/b**5/rho**(10.d0/3.d0))*hx**2
      vx=c*(v1*jj+v2*jja+v3*jjaa)
      end
