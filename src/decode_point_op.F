      subroutine decode_point_op(t,nt,i,a,djgen,lmax,jmax,su1,su1_j,fac)
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      character*80, intent(in) :: t
      integer, intent(in) :: nt
      integer :: i,k,l,m,m1,lmax,jmax,nrot,iwig,ist,ii,iwig0,jj,im,im1,
     &           ik
      real*8 :: v(3),a(maxwig),su1(-lmax:lmax,-lmax:lmax),fac(0:*),sp,
     &          v12,teta,fi,om,ddot,su,ss,wign
      complex*16 :: wign_j,djgen(maxwig_j),su1_j(jmax+1,jmax+1),cu,
     &              suj
      if(t(i:i).eq.'R') then
        i=i+1
        read(t(i:i),*)nrot
        i=i+1
        call parsvc(t,nt,i,v,rbas)
        sp=ddot(3,v,1,v,1)
        sp=1.d0/dsqrt(sp)
        do k=1,3
          v(k)=v(k)*sp
        enddo
        v12=v(1)**2+v(2)**2
        if(v12.gt.1.d-12) then
          teta=acos(v(3))
          fi=acos(v(1)/sqrt(v12))
          if(v(2).lt.0.d0) fi=-fi
        else
          teta=0.d0
          fi=0.d0
        endif
        om=2.d0*pi/nrot
        do l=0,lmax
          do m=-l,l
            do m1=-l,l
              su1(m,m1)=wign(l,m,m1,fi,teta,0.d0,fac)
            enddo
          enddo
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              su=0.d0
              do k=-l,l
                ss=k*om
                su=su+su1(m1,k)*(su1(m,k)*dcos(ss)+su1(m,-k)*dsin(ss))
              enddo
              a(iwig)=su
            enddo
          enddo
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
              jj=l+l+ii
              do m=-jj,jj,2
                im=(m+jj)/2+1
                do m1=-jj,jj,2
                  im1=(m1+jj)/2+1
                  su1_j(im,im1)=wign_j(jj,m,m1,fi,teta,0.d0,fac)
                enddo
              enddo
              do m=-jj,jj,2
                im=(m+jj)/2+1
                do m1=-jj,jj,2
                  im1=(m1+jj)/2+1
                  iwig=iwig0+(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
                  suj=(0.d0,0.d0)
                  do k=-jj,jj,2
                    ik=(k+jj)/2+1
                    ss=k*om/2.d0
                    cu=dcmplx(0.d0,ss)
                    cu=exp(cu)
                    suj=suj+dconjg(su1_j(im1,ik))*su1_j(im,ik)*cu
                  enddo
                  djgen(iwig)=suj
                enddo
              enddo
            enddo
          enddo
        endif
      else if(t(i:i).eq.'M') then
        i=i+1
        call parsvc(t,nt,i,v,rbas)
        sp=ddot(3,v,1,v,1)
        sp=1.d0/dsqrt(sp)
        do k=1,3
          v(k)=v(k)*sp
        enddo
        v12=v(1)**2+v(2)**2
        if(v12.gt.1.d-12) then
          teta=acos(v(3))
          fi=acos(v(1)/sqrt(v12))
          if(v(2).lt.0.d0) fi=-fi
        else
          teta=0.d0
          fi=0.d0
        endif
        do l=0,lmax
          do m=-l,l
            do m1=-l,l
              su1(m,m1)=wign(l,m,m1,fi,teta,0.d0,fac)
            enddo
          enddo
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              su=0.d0
              do k=-l,l
                su=su+(-1.d0)**(l+k)*su1(m1,k)*su1(m,k)
              enddo
              a(iwig)=su
            enddo
          enddo
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
              jj=l+l+ii
              do m=-jj,jj,2
                im=(m+jj)/2+1
                do m1=-jj,jj,2
                  im1=(m1+jj)/2+1
                  su1_j(im,im1)=wign_j(jj,m,m1,fi,teta,0.d0,fac)
                enddo
              enddo
              do m=-jj,jj,2
                im=(m+jj)/2+1
                do m1=-jj,jj,2
                  im1=(m1+jj)/2+1
                  iwig=iwig0+(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
                  suj=(0.d0,0.d0)
                  do k=-jj,jj,2
                    ik=(k+jj)/2+1
                    suj=suj+dconjg(su1_j(im1,ik))*su1_j(im,ik)
     &                         *(-1.d0)**(l+(k+1)/2)*(0.d0,1.d0)
                  enddo
                  djgen(iwig)=suj
                enddo
              enddo
            enddo
          enddo
        endif
      else if(t(i:i).eq.'I') then
        do l=0,lmax
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              if (m.eq.m1)then
                a(iwig)=(-1)**l
              else
                a(iwig)=0.d0
              endif
            enddo
          enddo
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
              jj=l+l+ii
              do m=-jj,jj,2
                do m1=-jj,jj,2
                  iwig=iwig0+(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
                  if (m.eq.m1)then
                    ss=(-1)**l
                  else
                    ss=0.d0
                  endif
                  djgen(iwig)=dcmplx(ss,0.d0)
                enddo
              enddo
            enddo
          enddo
        endif
        i=i+1
      else if(t(i:i).eq.'E') then
        do l=0,lmax
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              if (m.eq.m1)then
                a(iwig)=1.d0
              else
                a(iwig)=0.d0
              endif
            enddo
          enddo
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
              jj=l+l+ii
              do m=-jj,jj,2
                do m1=-jj,jj,2
                  iwig=iwig0+(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
                  if (m.eq.m1)then
                    ss=1.d0
                  else
                    ss=0.d0
                  endif
                  djgen(iwig)=dcmplx(ss,0.d0)
                enddo
              enddo
            enddo
          enddo
        endif
        i=i+1
      else
        if(maswrk) write(iun,*)'*** op must start with "R","M" or "I"'
        call ending
      endif
      end
