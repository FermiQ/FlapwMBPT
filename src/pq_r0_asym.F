      subroutine pq_r0_asym(m,l,ii,r,ai,bi,p,q,z)
      use manager_mod
      use units_mod
      implicit none
      integer, intent(in) :: m,l,ii
      real*8, intent(in) :: r(0:m),ai(0:3),bi(0:3),z
      real*8, intent(out) :: p(0:m),q(0:m)
      integer :: i,kap
      real*8 :: rll,lam,c2
C   ASSYMPTOTIC R ----> 0   BEHAVIOR ::
      c2=clight**2
      if(irel==0) then
        do i=0,m
          rll=r(i)**(l+1)
          p(i)=rll*(ai(0)+r(i)*(ai(1)+r(i)*(ai(2)+r(i)*ai(3))))
          q(i)=rll*(bi(0)+r(i)*(bi(1)+r(i)*(bi(2)+r(i)*bi(3))))
        enddo
      else if(irel==1) then
        if(z<1.d-3) then
          if(l==0) lam=0.d0
          if(l/=0) lam=-0.5d0+sqrt(1.d0+4*l*(l+1))/2
        else
          lam=sqrt(1.d0+l*(l+1)-4*z*z/c2)
        endif
        do i=0,m
          rll=r(i)**lam
          p(i)=rll*(ai(0)+ai(1)*r(i)+ai(2)*r(i)**2+ai(3)*r(i)**3)
          q(i)=rll*(bi(0)+bi(1)*r(i)+bi(2)*r(i)**2+bi(3)*r(i)**3)
        enddo
      else if(irel==2) then
        if(ii==-1) kap=l
        if(ii==1) kap=-l-1
        if(z>1.d-3) then
          lam=sqrt(kap*kap-4*z*z/c2)
        else
          if(ii==-1) lam=l
          if(ii==1) lam=l+1
        endif
        do i=0,m
          rll=r(i)**lam
          p(i)=rll*(ai(0)+r(i)*(ai(1)+r(i)*(ai(2)+r(i)*ai(3))))
          q(i)=rll*(bi(0)+r(i)*(bi(1)+r(i)*(bi(2)+r(i)*bi(3))))
        enddo
      endif
      end
      
      
      
      
      subroutine dpq_r0_asym(m,l,ii,r,ai,bi,p,q,z)
      use manager_mod
      use units_mod
      implicit none
      integer, intent(in) :: m,l,ii
      real*8, intent(in) :: r(0:m),ai(0:3),bi(0:3),z
      real*8, intent(out) :: p(0:m),q(0:m)
      integer :: i,kap
      real*8 :: rll,f,g,lam,c2
C   ASSYMPTOTIC R ----> 0   BEHAVIOR ::
      c2=clight**2
      if(irel==0) then
        do i=0,m
          f=ai(1)+2*ai(2)*r(i)+3*ai(3)*r(i)**2
          g=bi(1)+2*bi(2)*r(i)+3*bi(3)*r(i)**2
          rll=r(i)**(l+1)
          p(i)=rll*f
          q(i)=rll*g
          f=ai(0)+ai(1)*r(i)+ai(2)*r(i)**2+ai(3)*r(i)**3
          g=bi(0)+bi(1)*r(i)+bi(2)*r(i)**2+bi(3)*r(i)**3
          rll=(l+1)*r(i)**l
          p(i)=p(i)+rll*f
          q(i)=q(i)+rll*g
        enddo
      else if(irel==1) then
        if(z<1.d-3) then
          if(l==0) lam=0.d0
          if(l/=0) lam=-0.5d0+sqrt(1.d0+4*l*(l+1))/2
        else
          lam=sqrt(1.d0+l*(l+1)-4*z*z/c2)
        endif
        do i=0,m
          if(i==0.and.l==0.and.z<1.d-3) cycle  
          rll=r(i)**lam
          f=ai(0)+ai(1)*r(i)+ai(2)*r(i)**2+ai(3)*r(i)**3
          g=ai(1)+2*ai(2)*r(i)+3*ai(3)*r(i)**2
          p(i)=rll*g
          if(i/=0.or.lam>=1.d0) p(i)=p(i)+lam*r(i)**(lam-1)*f
          f=bi(0)+bi(1)*r(i)+bi(2)*r(i)**2+bi(3)*r(i)**3
          g=bi(1)+2*bi(2)*r(i)+3*bi(3)*r(i)**2
          q(i)=rll*g
          if(i/=0.or.lam>=1.d0) q(i)=q(i)+lam*r(i)**(lam-1)*f
        enddo
      else if(irel==2) then
        if(ii==-1) kap=l
        if(ii==1) kap=-l-1
        if(z>1.d-3) then
          lam=sqrt(kap*kap-4*z*z/c2)
        else
          if(ii==-1) lam=l
          if(ii==1) lam=l+1
        endif
        do i=0,m
          f=ai(1)+2*ai(2)*r(i)+3*ai(3)*r(i)**2
          g=bi(1)+2*bi(2)*r(i)+3*bi(3)*r(i)**2
          rll=r(i)**lam
          p(i)=rll*f
          q(i)=rll*g
          if(i/=0) then
            f=ai(0)+ai(1)*r(i)+ai(2)*r(i)**2+ai(3)*r(i)**3
            g=bi(0)+bi(1)*r(i)+bi(2)*r(i)**2+bi(3)*r(i)**3
            rll=lam*r(i)**(lam-1)
            p(i)=p(i)+rll*f
            q(i)=q(i)+rll*g
          endif
        enddo
      endif
      end