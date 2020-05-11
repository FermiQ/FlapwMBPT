      subroutine tensor_prod_via_ylm(t,t1,tf,lmax)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmax
      complex*16, intent(in) :: t(3,3),t1(3,3)
      complex*16, intent(out) :: tf(3,3)
      integer :: lmmax,lm,lm1,lm2,ii,iii,indx,icg1,icg2,icg
      complex*16 :: ty(9),t1y(9)
      complex*16, allocatable :: ry(:)
      lmmax=(lmax+1)**2
      call tensor_y_exp(t,ty)
      call tensor_y_exp(t1,t1y)
      allocate(ry(lmmax))
      ry=(0.d0,0.d0)
      do lm=1,9
        do lm1=1,9
          ii = max0(lm,lm1)
          iii = min0(lm,lm1)
          indx = (ii*(ii-1))/2 + iii
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do icg = icg1, icg2
            lm2 = jcg(icg)
            if(lm2>lmmax) cycle
            ry(lm2)=ry(lm2)+ty(lm)*t1y(lm1)*cg(icg)
          enddo
        enddo
      enddo
      call tensor_from_ylm(lmax,tf,ry)
      deallocate(ry)
      end


      subroutine tensor_prod_via_ylm_r(t,t1,tf,lmax)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmax
      real*8, intent(in) :: t(3,3),t1(3,3)
      real*8, intent(out) :: tf(3,3)
      integer :: lmmax,lm,lm1,lm2,ii,iii,indx,icg1,icg2,icg
      real*8 :: ty(9),t1y(9)
      real*8, allocatable :: ry(:)
      lmmax=(lmax+1)**2
      call tensor_y_exp_r(t,ty)
      call tensor_y_exp_r(t1,t1y)
      allocate(ry(lmmax))
      ry=0.d0
      do lm=1,9
        do lm1=1,9
          ii = max0(lm,lm1)
          iii = min0(lm,lm1)
          indx = (ii*(ii-1))/2 + iii
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do icg = icg1, icg2
            lm2 = jcg(icg)
            if(lm2>lmmax) cycle
            ry(lm2)=ry(lm2)+ty(lm)*t1y(lm1)*cg(icg)
          enddo
        enddo
      enddo
      call tensor_from_ylm_r(lmax,tf,ry)
      deallocate(ry)
      end
