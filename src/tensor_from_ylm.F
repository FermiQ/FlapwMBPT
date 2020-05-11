      subroutine tensor_from_ylm(lmax,t,ty)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmax
      complex*16, intent(in) :: ty((lmax+1)**2)
      complex*16, intent(out) :: t(3,3)
      integer :: it,lm
      real*8 :: sq2,dir(3)
      complex*16 :: t6(6)
      real*8, allocatable :: ylm(:)
      sq2=sqrt(2.d0)/2.d0
      allocate(ylm((lmax+1)**2))
      t6=(0.d0,0.d0)
      do it=1,6
        if(it==1) dir=(/1.d0,0.d0,0.d0/)
        if(it==2) dir=(/0.d0,1.d0,0.d0/)
        if(it==3) dir=(/0.d0,0.d0,1.d0/)
        if(it==4) dir=(/sq2,sq2,0.d0/)
        if(it==5) dir=(/sq2,0.d0,sq2/)
        if(it==6) dir=(/0.d0,sq2,sq2/)
        call sphharm(dir(1),dir(2),dir(3),lmax,ylm)
        do lm=1,(lmax+1)**2
          t6(it)=t6(it)+ty(lm)*ylm(lm)
        enddo
      enddo
      deallocate(ylm)
      call tensor_unpack(t6,t)
      end


      subroutine tensor_from_ylm_r(lmax,t,ty)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmax
      real*8, intent(in) :: ty((lmax+1)**2)
      real*8, intent(out) :: t(3,3)
      integer :: it,lm
      real*8 :: sq2,dir(3)
      complex*16 :: t6(6)
      real*8, allocatable :: ylm(:)
      sq2=sqrt(2.d0)/2.d0
      allocate(ylm((lmax+1)**2))
      t6=0.d0
      do it=1,6
        if(it==1) dir=(/1.d0,0.d0,0.d0/)
        if(it==2) dir=(/0.d0,1.d0,0.d0/)
        if(it==3) dir=(/0.d0,0.d0,1.d0/)
        if(it==4) dir=(/sq2,sq2,0.d0/)
        if(it==5) dir=(/sq2,0.d0,sq2/)
        if(it==6) dir=(/0.d0,sq2,sq2/)
        call sphharm(dir(1),dir(2),dir(3),lmax,ylm)
        do lm=1,(lmax+1)**2
          t6(it)=t6(it)+ty(lm)*ylm(lm)
        enddo
      enddo
      deallocate(ylm)
      call tensor_unpack_r(t6,t)
      end
