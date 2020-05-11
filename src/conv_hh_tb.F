      real*8 function conv_hh_tb(a,b,n)
c     calculates Sum_ij  a_ij*b_ji for two hermition matrices
c     They are stored in upper triangle packed form
      use parallel_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: a(n,n)
      complex*16, intent(in) :: b(nd_tb_ij(me_t_b+1))
      integer :: i,k0,k1,ij,j,ij0
      real*8 :: d
      k0=n0_tb_ij(me_t_b+1)
      k1=k0+nd_tb_ij(me_t_b+1)
      d=0.d0
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij<=k0) cycle
          if(ij>k1) cycle
          ij0=ij-k0
          if(i==j) d=d+a(i,i)*b(ij0)
          if(i/=j) d=d+2.d0*real(dcmplx(a(i,j),-a(j,i))*b(ij0))
        enddo
      enddo
      if(nproc_t*nproc_b/=1) call dgop(d,1,'  +',comm_t_b)
      conv_hh_tb=d
      end
