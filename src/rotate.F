      subroutine rotate(x2,y2,z2,x1,y1,z1,u,key)
	  use units_mod
      implicit none
	  integer :: key
      real*8 :: u(3,3),yg(3),yg1(3),r,x1,y1,z1,x2,y2,z2,s3
      r=x2*x2+y2*y2+z2*z2
      if(r<1.d-13)then
        x1=0.d0
        y1=0.d0
        z1=0.d0
        return
      endif
      s3=sqrt(0.75d0/pi/r)
      yg(1)=s3*y2
      yg(2)=s3*z2
      yg(3)=s3*x2
      if(key.eq.1)call prd1(yg,u,yg1) !!! inversional ig-operation
      if(key.eq.2)call prd2(yg,u,yg1) !!!      direct ig-operation
      x1=yg1(3)/s3
      y1=yg1(1)/s3
      z1=yg1(2)/s3
      end
      
      
      
	  subroutine rotate1(x,y,l,u0,key)
	  use solid_mod
      implicit none
c	  key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c	  key=2 Y(m)=Sum_m1 D_m1_m X_m1
	  integer :: l,key,nls,lm,m,nls1,m1,l1m1,iwig,l1,l2m,ll2m
      real*8 x(2*l+1),y(2*l+1),u0(maxwig),sum
	  l1=l+1
	  l2m=l1+l
	  ll2m=l2m*l
      nls=(l2m-2)*ll2m/3+l1+ll2m
	  lm=0
	  if(key.eq.1) then
        do m=-l,l
          lm=lm+1
          nls1=nls+m
	      sum=0.d0
          do m1=-l,l
            l1m1=l1+m1
            iwig=nls1+l2m*m1
            sum=sum+u0(iwig)*x(l1m1)
	      enddo
          y(lm)=sum
	    enddo
	  else if(key.eq.2) then
        do m=-l,l
          lm=lm+1
          nls1=nls+l2m*m
	      sum=0.d0
          do m1=-l,l
            l1m1=l1+m1
            iwig=nls1+m1
            sum=sum+u0(iwig)*x(l1m1)
	      enddo
          y(lm)=sum
	    enddo
	  endif
      end
      
      
	  subroutine rotate1c(x,y,l,u0,key)
	  use solid_mod
      implicit none
c	  key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c	  key=2 Y(m)=Sum_m1 D_m1_m X_m1
	  integer :: key,l,nls,lm,m,nls1,m1,l1m1,iwig
      real*8 u0(maxwig)
      complex*16 x(2*l+1),y(2*l+1),sum
      nls=l*(2*l-1)*(2*l+1)/3+l+1
	  lm=0
	  if(key.eq.1) then
        do m=-l,l
          lm=lm+1
          nls1=nls+m
	      sum=(0.d0,0.d0)
          do m1=-l,l
            l1m1=l+m1+1
            iwig=nls1+(2*l+1)*(l+m1)
            sum=sum+u0(iwig)*x(l1m1)
	      enddo
          y(lm)=sum
	    enddo
	  else if(key.eq.2) then
        do m=-l,l
          lm=lm+1
          nls1=nls+(2*l+1)*(l+m)
	      sum=(0.d0,0.d0)
          do m1=-l,l
            l1m1=l+m1+1
            iwig=nls1+m1
            sum=sum+u0(iwig)*x(l1m1)
	      enddo
          y(lm)=sum
	    enddo
	  endif
      end
      
      
      
	  subroutine rotate2(x,y,l,ii,u0,key,key1)
	  use solid_mod
      implicit none
c	  key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c	  key=2 Y(m)=Sum_m1 D_m1_m X_m1
c	  key1=1 conjugated rotation  (D*)
c	  key1=2 usual                (D )
	  integer :: key,key1,l,ii,jj,iwig0,lm,m,nls1,l1m1,m1,iwig,i5
      complex*16 x(*),y(*),u0(maxwig_j),sum
	  jj=l+l+ii
	  if(key1.eq.1) then
	    iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
        i5=l+(ii+1)/2
	    if(key.eq.1) then
	      lm=0
          do m=-jj,jj,2
	        lm=lm+1
            nls1=iwig0+(ii+m)/2+i5*ii
	        sum=(0.d0,0.d0)
	        l1m1=0
            do m1=-jj,jj,2
	          l1m1=l1m1+1
              iwig=nls1+i5*m1
              sum=sum+conjg(u0(iwig))*x(l1m1)
	        enddo
            y(lm)=sum
	      enddo
	    else if(key.eq.2) then
	      lm=0
          do m=-jj,jj,2
	        lm=lm+1
            nls1=iwig0+(2*l+ii+1)*(ii+m)/2
	        sum=(0.d0,0.d0)
	        l1m1=0
            do m1=-jj,jj,2
	          l1m1=l1m1+1
              iwig=nls1+(ii+m1)/2
              sum=sum+conjg(u0(iwig))*x(l1m1)
	        enddo
            y(lm)=sum
	      enddo
	    endif
	  else if(key1.eq.2) then
	    iwig0=(8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
	    if(key.eq.1) then
	      lm=0
          do m=-jj,jj,2
	        lm=lm+1
            nls1=iwig0+(ii+m)/2
	        sum=(0.d0,0.d0)
	        l1m1=0
            do m1=-jj,jj,2
	          l1m1=l1m1+1
              iwig=nls1+(2*l+ii+1)*(ii+m1)/2
              sum=sum+u0(iwig)*x(l1m1)
	        enddo
            y(lm)=sum
	      enddo
	    else if(key.eq.2) then
	      lm=0
          do m=-jj,jj,2
	        lm=lm+1
            nls1=iwig0+(2*l+ii+1)*(ii+m)/2
	        sum=(0.d0,0.d0)
	        l1m1=0
            do m1=-jj,jj,2
	          l1m1=l1m1+1
              iwig=nls1+(ii+m1)/2
              sum=sum+u0(iwig)*x(l1m1)
	        enddo
            y(lm)=sum
	      enddo
	    endif
	  endif
      end
      
      
      subroutine rotate2d(x2,y2,x1,y1,u,key)
	  use units_mod
      implicit none
	  integer :: key
      real*8 :: u(2,2),x1,y1,x2,y2
      if(key==1) then
        x1=x2*u(1,1)+y2*u(2,1)
        y1=x2*u(1,2)+y2*u(2,2)
      else if(key==2) then
        x1=u(1,1)*x2+u(1,2)*y2
        y1=u(2,1)*x2+u(2,2)*y2
      endif
      end
      
      
      subroutine rotate_c(v0,v,u,key)
	  use units_mod
      implicit none
	  integer, intent(in) :: key
      real*8, intent(in) :: u(3,3)
      complex*16, intent(in) :: v0(3)
      complex*16, intent(out) :: v(3)
      integer :: i
      real*8 :: yg(3),yg1(3),r,s3,x0(3),x(3)
      do i=1,2
        if(i==1) x0=real(v0)
        if(i==2) x0=imag(v0)
        r=dot_product(x0,x0)
        if(r==0.d0) then
          x=0.d0
        else
          s3=sqrt(0.75d0/pi/r)
          yg(1)=s3*x0(2)
          yg(2)=s3*x0(3)
          yg(3)=s3*x0(1)
          if(key.eq.1)call prd1(yg,u,yg1) !!! inversional ig-operation
          if(key.eq.2)call prd2(yg,u,yg1) !!!      direct ig-operation
          x(1)=yg1(3)/s3
          x(2)=yg1(1)/s3
          x(3)=yg1(2)/s3
        endif
        if(i==1) v=x
        if(i==2) v=v+(0.d0,1.d0)*x
      enddo
      end
      
      
      subroutine rotate_r(v0,v,u,key)
	  use units_mod
      implicit none
	  integer, intent(in) :: key
      real*8, intent(in) :: u(3,3),v0(3)
      real*8, intent(out) :: v(3)
      real*8 :: yg(3),yg1(3),r,s3
      r=dot_product(v0,v0)
      if(r==0.d0) then
        v=0.d0
      else
        s3=sqrt(0.75d0/pi/r)
        yg(1)=s3*v0(2)
        yg(2)=s3*v0(3)
        yg(3)=s3*v0(1)
        if(key.eq.1)call prd1(yg,u,yg1) !!! inversional ig-operation
        if(key.eq.2)call prd2(yg,u,yg1) !!!      direct ig-operation
        v(1)=yg1(3)/s3
        v(2)=yg1(1)/s3
        v(3)=yg1(2)/s3
      endif
      end
