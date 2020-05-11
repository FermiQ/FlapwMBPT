      subroutine pack_fxc(n,f,fp)
      use parallel_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      complex*16, intent(in) :: f(n,n)
      complex*16, intent(out) :: fp(nd_tb_fxc(me_t_b+1))
      integer :: i,j,k0,k1,ij
      k0=n0_tb_fxc(me_t_b+1)
      k1=k0+nd_tb_fxc(me_t_b+1)
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij<=k0) cycle
          if(ij>k1) cycle
          fp(ij-k0)=f(i,j)
        enddo
      enddo
      end


      subroutine pack_fxc_r(n,f,fp)
      use parallel_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: f(n,n)
      real*8, intent(out) :: fp(nd_tb_fxc(me_t_b+1))
      integer :: i,j,k0,k1,ij
      k0=n0_tb_fxc(me_t_b+1)
      k1=k0+nd_tb_fxc(me_t_b+1)
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij<=k0) cycle
          if(ij>k1) cycle
          fp(ij-k0)=f(i,j)
        enddo
      enddo
      end


      subroutine unpack_fxc(n,f,fp)
      use parallel_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      complex*16, intent(in) :: fp(nd_tb_fxc(me_t_b+1))
      complex*16, intent(out) :: f(n,n)
      integer :: i,j,k0,k1,ij
      k0=n0_tb_fxc(me_t_b+1)
      k1=k0+nd_tb_fxc(me_t_b+1)
      f=(0.d0,0.d0)
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij<=k0) cycle
          if(ij>k1) cycle
          f(i,j)=fp(ij-k0)
          if(i/=j) f(j,i)=conjg(fp(ij-k0))
        enddo
      enddo
      if(nproc_t*nproc_b/=1) call dgop(f,2*n*n,'  +',comm_t_b)
      end


      subroutine unpack_fxc_r(n,f,fp)
      use parallel_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: fp(nd_tb_fxc(me_t_b+1))
      real*8, intent(out) :: f(n,n)
      integer :: i,j,k0,k1,ij
      k0=n0_tb_fxc(me_t_b+1)
      k1=k0+nd_tb_fxc(me_t_b+1)
      f=0.d0
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij<=k0) cycle
          if(ij>k1) cycle
          f(i,j)=fp(ij-k0)
          if(i/=j) f(j,i)=fp(ij-k0)
        enddo
      enddo
      if(nproc_t*nproc_b/=1) call dgop(f,n*n,'  +',comm_t_b)
      end
