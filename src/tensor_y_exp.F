      subroutine tensor_y_exp(t,ty)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: t(3,3)
      complex*16, intent(out) :: ty(9)
      integer :: i,j,lm,lm1,lm2,ii,iii,indx,icg1,icg2,icg
      ty=(0.d0,0.d0)
      do j=1,3
        if(j==1) lm1=2
        if(j==2) lm1=4
        if(j==3) lm1=3
        do i=1,3
          if(i==1) lm=2
          if(i==2) lm=4
          if(i==3) lm=3
          ii = max0(lm,lm1)
          iii = min0(lm,lm1)
          indx = (ii*(ii-1))/2 + iii
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do icg = icg1, icg2
            lm2 = jcg(icg)
            ty(lm2)=ty(lm2)+t(i,j)*cg(icg)
          enddo
        enddo
      enddo
      ty=4.d0*pi/3.d0*ty
      end


      subroutine tensor_y_exp_r(t,ty)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) :: t(3,3)
      real*8, intent(out) :: ty(9)
      integer :: i,j,lm,lm1,lm2,ii,iii,indx,icg1,icg2,icg
      ty=0.d0
      do j=1,3
        if(j==1) lm1=2
        if(j==2) lm1=4
        if(j==3) lm1=3
        do i=1,3
          if(i==1) lm=2
          if(i==2) lm=4
          if(i==3) lm=3
          ii = max0(lm,lm1)
          iii = min0(lm,lm1)
          indx = (ii*(ii-1))/2 + iii
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do icg = icg1, icg2
            lm2 = jcg(icg)
            ty(lm2)=ty(lm2)+t(i,j)*cg(icg)
          enddo
        enddo
      enddo
      ty=4.d0*pi/3.d0*ty
      end
