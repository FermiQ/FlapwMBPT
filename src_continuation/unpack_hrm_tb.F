      subroutine unpack_hrm_tb(aunpck,apck,n,n0)
      use parallel_mod
      implicit none
      integer, intent(in) :: n,n0
      complex*16, intent(in) :: apck(nd_tb_ij(me_t_b+1))
      complex*16, intent(out) :: aunpck(n0,n0)
      integer :: i,j,k0,k1,ij
      k0=n0_tb_ij(me_t_b+1)
      k1=k0+nd_tb_ij(me_t_b+1)
      aunpck=(0.d0,0.d0)
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij<=k0) cycle
          if(ij>k1) cycle
          aunpck(i,j)=apck(ij-k0)
          if(i/=j) aunpck(j,i)=conjg(aunpck(i,j))
        enddo
      enddo
      if(nproc_t*nproc_b/=1) call dgop(aunpck,2*n0*n0,'  +',comm_t_b)
      end
