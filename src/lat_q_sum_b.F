      subroutine lat_q_sum_b(q,idt,lmax,dl)
c     k-space part of reduced structure constants (Ewald).
c     tau in alat-units. q and rlat in tpiba-units.
c ------- We calculate only regular contribution here --------------
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: idt,lmax
      complex*16, intent(out) :: dl((lmax+1)**2)
      integer :: nlm,ilm,ir,l,m
      real*8 :: q(3),dta(3),r(3),tpi,gamma,fpibv,tpiba,qq,scalp,r2,r1,
     &          den0,sqpi4
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
      allocate(yl(nlm))
      do ir=1,nkr
        r=tpiba*(q+rlat(:,ir))
        r2=dot_product(r,r)
        if(r2<1.d-12) then
          dl(1)=dl(1)+exp(-gamma*omega_hse**2)/omega_hse**2/sqpi4
        else
          scalp=par*dot_product(r,dta)
          eiphi=dcmplx(dcos(scalp),dsin(scalp))
          r1=dsqrt(r2)
          call sphharm(r(1)/r1,r(2)/r1,r(3)/r1,lmax,yl)
          r2=r2+omega_hse**2
          den0=exp(-gamma*r2)/r2
          ilm=0
          do l=0,lmax
            cs=eiphi*den0*r1**l
            do m=-l,l
              ilm=ilm+1
              dl(ilm)=dl(ilm)+yl(ilm)*cs
            enddo
          enddo
        endif
      enddo
      eiphi=(0.d0,1.d0)*fpibv
      ilm=0
      do l=0,lmax
        eiphi=eiphi*(0.d0,-1.d0)
        do m=-l,l
          ilm=ilm+1
          dl(ilm)=eiphi*dl(ilm)
        enddo
      enddo
      deallocate(yl)
      end
