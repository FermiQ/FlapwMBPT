      subroutine vec_tens_prod(t,f)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: t(3,3)
      complex*16, intent(out) :: f(3,3)
      integer :: lm,lm1,lm2,m,m1,ii,iii,indx,icg1,icg2,icg,i,j
      complex*16, allocatable :: ty(:)
      allocate(ty(9))
      f=(0.d0,0.d0)
      call tensor_y_exp(t,ty)
      do i=1,3
        if(i==1) m=1
        if(i==2) m=-1
        if(i==3) m=0
        lm=3+m
        do j=1,3
          if(j==1) m1=1
          if(j==2) m1=-1
          if(j==3) m1=0
          lm1=3+m1
          ii = max0(lm,lm1)
          iii = min0(lm,lm1)
          indx = (ii*(ii-1))/2 + iii
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do icg = icg1, icg2
            lm2 = jcg(icg)
            f(i,j)=f(i,j)+ty(lm2)*cg(icg)
          enddo
        enddo
      enddo
      deallocate(ty)
      end

      subroutine vec_tens_prod_r(t,f)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) :: t(3,3)
      real*8, intent(out) :: f(3,3)
      integer :: lm,lm1,lm2,m,m1,ii,iii,indx,icg1,icg2,icg,i,j
      real*8, allocatable :: ty(:)
      allocate(ty(9))
      f=0.d0
      call tensor_y_exp_r(t,ty)
      do i=1,3
        if(i==1) m=1
        if(i==2) m=-1
        if(i==3) m=0
        lm=3+m
        do j=1,3
          if(j==1) m1=1
          if(j==2) m1=-1
          if(j==3) m1=0
          lm1=3+m1
          ii = max0(lm,lm1)
          iii = min0(lm,lm1)
          indx = (ii*(ii-1))/2 + iii
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do icg = icg1, icg2
            lm2 = jcg(icg)
            f(i,j)=f(i,j)+ty(lm2)*cg(icg)
          enddo
        enddo
      enddo
      deallocate(ty)
      end
