      subroutine spline_matrix(n,hi,key,aa,bb,cc,dd)
c     key = 0  - Free boundary conditions
c     key = 2  - First order derivatives at boundaries are given      
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,key
      real*8, intent(in) :: hi(n)
      real*8, intent(out) :: aa(n,0:n+key),bb(n,0:n+key),cc(n,0:n+key),
     &                       dd(n,0:n+key)
      integer :: i,j
      real*8, allocatable :: a0(:,:),r0(:,:)
c ---------------------------------------------------------------
      if(key==0) then
        aa=0.d0
        do i=1,n
          aa(i,i-1)=1.d0
        enddo
        allocate(a0(2:n,2:n))
        allocate(r0(2:n,0:n))
        a0=0.d0
        r0=0.d0
        do i=2,n
          if(i>2) a0(i,i-1)=-hi(i-1)/3
          a0(i,i)=-2.d0*(hi(i-1)+hi(i))/3.d0
          if(i<n) a0(i,i+1)=-hi(i)/3
          r0(i,i-2)=-1.d0/hi(i-1)
          r0(i,i-1)=1.d0/hi(i-1)+1.d0/hi(i)
          r0(i,i)=-1.d0/hi(i)
        enddo
        call invers_r(n-1,a0,n-1)
        call dgemm('n','n',n-1,n+1,n-1,1.d0,a0,n-1,r0(2,0),n-1,0.d0,
     &             cc(2,0),n)
        cc(1,:)=0.d0
        deallocate(a0,r0)
c ---------------------------------------------------------------
        bb(1,:)=-hi(1)*cc(2,:)/3.d0
        bb(1,0)=bb(1,0)-1.d0/hi(1)
        bb(1,1)=bb(1,1)+1.d0/hi(1)
        do i=1,n-1
          bb(i+1,:)=bb(i,:)+hi(i)*(cc(i,:)+cc(i+1,:))
        enddo
c ---------------------------------------------------------------
        do i=1,n-1
          dd(i,:)=(cc(i+1,:)-cc(i,:))/3.d0/hi(i)
        enddo
        dd(n,:)=-cc(n,:)/3.d0/hi(n)
      else if(key==2) then
        aa=0.d0
        do i=1,n
          aa(i,i-1)=1.d0
        enddo
        allocate(a0(2:n,2:n))
        allocate(r0(2:n,0:n))
        a0=0.d0
        r0=0.d0
        a0(2,2)=-(hi(1)/2+2.d0*hi(2)/3.d0)
        a0(2,3)=-hi(2)/3
        r0(2,0)=-1.5d0/hi(1)
        r0(2,1)=1.5d0/hi(1)+1.d0/hi(2)
        r0(2,2)=-1.d0/hi(2)
        do i=3,n-1
          a0(i,i-1)=-hi(i-1)/3
          a0(i,i)=-2.d0*(hi(i-1)+hi(i))/3.d0
          a0(i,i+1)=-hi(i)/3
          r0(i,i-2)=-1.d0/hi(i-1)
          r0(i,i-1)=1.d0/hi(i-1)+1.d0/hi(i)
          r0(i,i)=-1.d0/hi(i)
        enddo
        a0(n,n-1)=-hi(n-1)/3
        a0(n,n)=-2.d0*hi(n-1)/3-hi(n)/2.d0
        r0(n,n-2)=-1.d0/hi(n-1)
        r0(n,n-1)=1.d0/hi(n-1)+1.5d0/hi(n)
        r0(n,n)=-1.5d0/hi(n)
        call invers_r(n-1,a0,n-1)
        call dgemm('n','n',n-1,n+1,n-1,1.d0,a0,n-1,r0(2,0),n-1,0.d0,
     &             cc(2,0),n)
        do i=2,n
          cc(i,n+1)=-0.5d0*a0(i,2)
          cc(i,n+2)=0.5d0*a0(i,n)
        enddo
        do j=0,n+2
          cc(1,j)=-0.5d0*cc(2,j)
        enddo
        cc(1,0)=-1.5d0/hi(1)**2+cc(1,0)
        cc(1,1)=1.5d0/hi(1)**2+cc(1,1)
        cc(1,n+1)=-1.5d0/hi(1)+cc(1,n+1)
        deallocate(a0,r0)
        bb(1,:)=0.d0
        bb(1,n+1)=1.d0
        do i=1,n-1
          bb(i+1,:)=bb(i,:)+hi(i)*(cc(i,:)+cc(i+1,:))
        enddo
        do i=1,n-1
          dd(i,:)=(cc(i+1,:)-cc(i,:))/3.d0/hi(i)
        enddo
        dd(n,:)=-(bb(n,:)+2.d0*hi(n)*cc(n,:))/3.d0/hi(n)**2
        dd(n,n+2)=dd(n,n+2)+1.d0/3.d0/hi(n)**2
      endif
      end