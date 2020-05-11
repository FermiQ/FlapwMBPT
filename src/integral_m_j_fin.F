      subroutine integral_m_j_fin(n,f,r1,iik,res)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: f(0:n,2),iik(n,0:3),r1(0:n)
      real*8, intent(out) :: res
      integer :: i
      real*8 :: a,f0,f1,df0,df1,a0,a1,a2,a3
      res=0.d0
      do i=1,n
        a=r1(i)-r1(i-1)
        f0=f(i-1,1)
        f1=f(i,1)
        df0=a*f(i-1,2)
        df1=a*f(i,2)
        a0=f0
        a1=df0
        a2=3*(f1-f0)-df1-2*df0
        a3=2*(f0-f1)+df0+df1
        res=res+a0*iik(i,0)+a1*iik(i,1)+a2*iik(i,2)+a3*iik(i,3)
      enddo
      end
