      subroutine tau_from_freq_spl(nt,tm,n,wn,tran,key)
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer, intent(in) :: nt,n,key
      real*8, intent(in) :: tm(0:nt),wn(0:n)
      real*8, intent(out) :: tran(0:n,0:nt,2)
      integer :: i,j,i_tau,m
      real*8 :: d,d2,d3,ct,w,pi2b,w0,a,b,xx
      real*8, allocatable :: aa(:,:),bb(:,:),cc(:,:),dd(:,:),
     &                       hi(:),ii(:,:,:),x(:)
      pi2b=(pi+pi)/betta_t
c ---------------------------------------------------------------
      allocate(x(0:n))
      if(key==0) then
        x=wn
      else if(key==1) then
        a=1.d0/wn(n)-1.d0/wn(0)
        b=1.d0/wn(0)
        do i=0,n
          x(i)=(1.d0/wn(i)-b)/a
        enddo
      else if(key==2) then
        a=1.d0/wn(n)**2-1.d0/wn(0)**2
        b=1.d0/wn(0)**2
        do i=0,n
          x(i)=(1.d0/wn(i)**2-b)/a
        enddo
      else if(key==3) then
        a=1.d0/wn(n)**3-1.d0/wn(0)**3
        b=1.d0/wn(0)**3
        do i=0,n
          x(i)=(1.d0/wn(i)**3-b)/a
        enddo
      endif
      allocate(hi(n))
      do i=1,n
        hi(i)=x(i)-x(i-1)
      enddo
      allocate(aa(n,0:n))
      allocate(bb(n,0:n))
      allocate(cc(n,0:n))
      allocate(dd(n,0:n))
      call spline_matrix(n,hi,0,aa,bb,cc,dd)
c ---------------------------------------------------------------
      allocate(ii(2:n-1,4,2))
      tran=0.d0
      do i_tau=0,nt
        ii=0.d0
        do i=2,n-1
          w0=wn(i-1)
          m=idnint((wn(i)-w0)/pi2b)
          do j=1,m
            d=pi2b*j
            w=w0+d
            if(key==1) then
              xx=(1.d0/w-b)/a
              d=xx-x(i-1)
            endif
            if(key==2) then
              xx=(1.d0/w/w-b)/a
              d=xx-x(i-1)
            endif
            if(key==3) then
              xx=(1.d0/w/w/w-b)/a
              d=xx-x(i-1)
            endif
            d2=d*d
            d3=d2*d
            ct=cos(w*tm(i_tau))
            ii(i,1,1)=ii(i,1,1)+ct
            ii(i,2,1)=ii(i,2,1)+d*ct
            ii(i,3,1)=ii(i,3,1)+d2*ct
            ii(i,4,1)=ii(i,4,1)+d3*ct
            ct=sin(w*tm(i_tau))
            ii(i,1,2)=ii(i,1,2)+ct
            ii(i,2,2)=ii(i,2,2)+d*ct
            ii(i,3,2)=ii(i,3,2)+d2*ct
            ii(i,4,2)=ii(i,4,2)+d3*ct
          enddo
        enddo
        do i=1,2
          do j=0,n
            tran(j,i_tau,i)=dot_product(aa(2:n-1,j),ii(:,1,i))
     &                     +dot_product(bb(2:n-1,j),ii(:,2,i))
     &                     +dot_product(cc(2:n-1,j),ii(:,3,i))
     &                     +dot_product(dd(2:n-1,j),ii(:,4,i))
          enddo
        enddo
      enddo   !! over i_tau
      deallocate(aa,bb,cc,dd,hi,x,ii)
      end
      
      
      
      
      
      subroutine dtau_from_freq_spl(n,wn,tran,key)
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,key
      real*8, intent(in) :: wn(0:n)
      real*8, intent(out) :: tran(0:n,2)
      integer :: i,j,m,mat
      real*8 :: d,d2,d3,ct,w,pi2b,w0,a,b,xx,pib
      real*8, allocatable :: aa(:,:),bb(:,:),cc(:,:),dd(:,:),
     &                       hi(:),ii(:,:,:),x(:)
      pib=pi/betta_t
      pi2b=pib+pib
c ---------------------------------------------------------------
      allocate(x(0:n))
      if(key==0) then
        x=wn
      else if(key==2) then
        a=1.d0/wn(n)**2-1.d0/wn(0)**2
        b=1.d0/wn(0)**2
        do i=0,n
          x(i)=(1.d0/wn(i)**2-b)/a
        enddo
      else if(key==3) then
        a=1.d0/wn(n)**3-1.d0/wn(0)**3
        b=1.d0/wn(0)**3
        do i=0,n
          x(i)=(1.d0/wn(i)**3-b)/a
        enddo
      endif
      allocate(hi(n))
      do i=1,n
        hi(i)=x(i)-x(i-1)
      enddo
      allocate(aa(n,0:n))
      allocate(bb(n,0:n))
      allocate(cc(n,0:n))
      allocate(dd(n,0:n))
      call spline_matrix(n,hi,0,aa,bb,cc,dd)
c ---------------------------------------------------------------
      allocate(ii(2:n-1,4,2))
      tran=0.d0
      ii=0.d0
      do i=2,n-1
        w0=wn(i-1)
        m=idnint((wn(i)-w0)/pi2b)
        do j=1,m
          d=pi2b*j
          w=w0+d
          mat=idnint((w-pib)/pi2b)
          if(key==2) then
            xx=(1.d0/w/w-b)/a
            d=xx-x(i-1)
          endif
          if(key==3) then
            xx=(1.d0/w/w/w-b)/a
            d=xx-x(i-1)
          endif
          d2=d*d
          d3=d2*d
          ct=-(-1)**mat*w
          ii(i,1,1)=ii(i,1,1)+ct
          ii(i,2,1)=ii(i,2,1)+d*ct
          ii(i,3,1)=ii(i,3,1)+d2*ct
          ii(i,4,1)=ii(i,4,1)+d3*ct
          ct=w
          ii(i,1,2)=ii(i,1,2)+ct
          ii(i,2,2)=ii(i,2,2)+d*ct
          ii(i,3,2)=ii(i,3,2)+d2*ct
          ii(i,4,2)=ii(i,4,2)+d3*ct
        enddo
      enddo
      do i=1,2
        do j=0,n
          tran(j,i)=dot_product(aa(2:n-1,j),ii(:,1,i))
     &             +dot_product(bb(2:n-1,j),ii(:,2,i))
     &             +dot_product(cc(2:n-1,j),ii(:,3,i))
     &             +dot_product(dd(2:n-1,j),ii(:,4,i))
        enddo
      enddo
      deallocate(aa,bb,cc,dd,hi,x,ii)
      end