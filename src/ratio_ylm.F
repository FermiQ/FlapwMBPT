      subroutine ratio_ylm(lmmax,t,ty)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: lmmax
      complex*16, intent(in) :: t(9)
      complex*16, intent(inout) :: ty(lmmax)
      integer :: lm,lm1,lm2,ii,iii,indx,icg1,icg2,info,icg
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: a(:,:)
      allocate(a(lmmax,lmmax))
      a=(0.d0,0.d0)
      do lm=1,lmmax
        do lm2=1,lmmax
          ii = max0(lm,lm2)
          iii = min0(lm,lm2)
          indx = (ii*(ii-1))/2 + iii
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do icg = icg1, icg2
            lm1 = jcg(icg)
            if(lm1>9) cycle
            a(lm2,lm)=a(lm2,lm)+t(lm1)*cg(icg)
          enddo
        enddo
      enddo
      allocate(ipiv(lmmax))
      call zgesv(lmmax,1,a,lmmax,ipiv,ty,lmmax,info)
      deallocate(a,ipiv)
      end
