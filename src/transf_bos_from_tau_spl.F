      subroutine transf_bos_from_tau_spl(nt,tm,nnu,wn,tran,n0,key)
c     key = 0  - Free boundary conditions
c     key = 2  - First order derivatives at boundaries are given
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer, intent(in) :: nt,nnu,key,n0
      real*8, intent(in) :: tm(0:nt),wn(0:nnu)
      real*8, intent(out) :: tran(0:nt+key,0:n0,2)
      integer :: i,j,i_nu
      real*8 :: c0,s0,c1,s1,c2,s2,c3,s3,w,wh,wt,cwt,swt,cwh,swh,wh2,wh3,
     &          wh4,wh5
      real*8, allocatable :: aa(:,:),bb(:,:),cc(:,:),dd(:,:),hi(:),
     &                       ii(:,:,:)
c ---------------------------------------------------------------
      allocate(hi(nt))
      do i=1,nt
        hi(i)=tm(i)-tm(i-1)
      enddo
      allocate(aa(nt,0:nt+key))
      allocate(bb(nt,0:nt+key))
      allocate(cc(nt,0:nt+key))
      allocate(dd(nt,0:nt+key))
      call spline_matrix(nt,hi,key,aa,bb,cc,dd)
c ---------------------------------------------------------------
      allocate(ii(nt,4,2))
      do i_nu=0,nnu
        w=wn(i_nu)
        do i=1,nt
          wh=w*hi(i)
          wt=w*tm(i-1)
          cwt=cos(wt)
          swt=sin(wt)
          cwh=cos(wh)
          swh=sin(wh)
          if(abs(w)<1.d-10) then
            c0=hi(i)
            c1=hi(i)**2/2
            c2=hi(i)**3/3
            c3=hi(i)**4/4
            s0=0.d0
            s1=0.d0
            s2=0.d0
            s3=0.d0
          else if(abs(wh)<1.2d-2) then
            wh2=wh*wh
            wh3=wh2*wh
            wh4=wh2*wh2
            wh5=wh3*wh2
            c0=hi(i)*(1.d0-wh2/6.d0+wh4/120.d0)
            c1=hi(i)**2*(0.5d0-wh2/8.d0+wh4/144.d0)
            c2=hi(i)**3*(1.d0/3.d0-wh2/10.d0+wh4/168.d0)
            c3=hi(i)**4*(0.25d0-wh2/12.d0+wh4/192.d0)
            s0=hi(i)*(wh/2.d0-wh3/24.d0+wh5/720.d0)
            s1=hi(i)**2*(wh/3.d0-wh3/30.d0+wh5/840.d0)
            s2=hi(i)**3*(wh/4.d0-wh3/36.d0+wh5/960.d0)
            s3=hi(i)**4*(wh/5.d0-wh3/42.d0+wh5/1080.d0)
          else
            c0=swh/w
            c1=(wh*swh+cwh-1.d0)/w**2
            c2=((wh*wh-2.d0)*swh+2.d0*wh*cwh)/w**3
            c3=(wh*(wh*wh-6.d0)*swh+3.d0*(wh*wh-2.d0)*cwh+6.d0)/w**4
            s0=(1.d0-cwh)/w
            s1=(swh-wh*cwh)/w**2
            s2=((2.d0-wh*wh)*cwh+2.d0*wh*swh-2.d0)/w**3
            s3=(wh*(6.d0-wh*wh)*cwh+3.d0*(wh*wh-2.d0)*swh)/w**4
          endif
          ii(i,1,1)=cwt*c0-swt*s0
          ii(i,2,1)=cwt*c1-swt*s1
          ii(i,3,1)=cwt*c2-swt*s2
          ii(i,4,1)=cwt*c3-swt*s3
          ii(i,1,2)=swt*c0+cwt*s0
          ii(i,2,2)=swt*c1+cwt*s1
          ii(i,3,2)=swt*c2+cwt*s2
          ii(i,4,2)=swt*c3+cwt*s3
        enddo
        do i=1,2
          do j=0,nt+key
            tran(j,i_nu,i)=dot_product(aa(:,j),ii(:,1,i))
     &                    +dot_product(bb(:,j),ii(:,2,i))
     &                    +dot_product(cc(:,j),ii(:,3,i))
     &                    +dot_product(dd(:,j),ii(:,4,i))
          enddo
        enddo
      enddo   !! over i_nu
      deallocate(aa,bb,cc,dd,ii)
      end