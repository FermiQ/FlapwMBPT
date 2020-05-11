      subroutine pack_hrm_tb(aunpck,apck,n,n0,key)
c     key = 0 creates new APCK
c     key = 1 adds to the already available APCK
      use parallel_mod
      implicit none
      integer, intent(in) :: n,n0,key
      complex*16, intent(in) :: aunpck(n0,n0)
      complex*16, intent(inout) :: apck(nd_tb_ij(me_t_b+1))
      integer :: i,j,k0,k1,ij
      k0=n0_tb_ij(me_t_b+1)
      k1=k0+nd_tb_ij(me_t_b+1)
      ij=0
      if(key==0) then
        do j=1,n
          do i=1,j
            ij=ij+1
            if(ij<=k0) cycle
            if(ij>k1) cycle
            apck(ij-k0)=aunpck(i,j)
          enddo
        enddo
      else if(key==1) then
        do j=1,n
          do i=1,j
            ij=ij+1
            if(ij<=k0) cycle
            if(ij>k1) cycle
            apck(ij-k0)=apck(ij-k0)+aunpck(i,j)
          enddo
        enddo
      endif
      end
