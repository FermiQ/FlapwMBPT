      subroutine phi_rpa_c_q0_met(n0,n,p,p00,p01)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0
      real*8, intent(in) :: p00
      complex*16, intent(in) :: p(n0,n0),p01(n)
      integer :: i,j,n1
      real*8 :: s,pi8,spi8
      real*8, allocatable :: d(:)
      complex*16, allocatable :: vp(:,:)
      n1=n+1
      pi8=8.d0*pi
      spi8=sqrt(pi8)
      allocate(vp(n1,n1))
      allocate(d(n1))
c --------------------- Trace of PV --------------------------------
      s=pi8*q2aver*p00
      do i=1,n
        s=s+p(i,i)*v_opt_e(i)
      enddo
c --------------------- Trace of LOG(1-VP) -------------------------
      s=s+log(pi8)-2.d0*logq_aver
      do i=1,n
        s=s+log(v_opt_e(i))
      enddo
      vp(1,1)=-p00
      vp(1,2:n1)=-p01
      vp(2:n1,1)=-conjg(p01)
      do j=2,n1
        do i=2,n1
          vp(i,j)=-p(i-1,j-1)
        enddo
        vp(j,j)=vp(j,j)+1.d0/v_opt_e(j-1)
      enddo
      call eig_val_solver(n1,n1,vp,d)
      do i=1,n1
        s=s+log(d(i))
      enddo
      deallocate(d,vp)
      s=s/2.d0
      phi_functional_c=phi_functional_c+et_sum_nu(0)*s*wgt(1)*2.d0
      end


      subroutine phi_rpa_c_q0_met_r(n0,n,p,p00,p01)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0
      real*8, intent(in) :: p00,p(n0,n0),p01(n)
      integer :: i,j,n1
      real*8 :: s,pi8,spi8
      real*8, allocatable :: d(:),vp(:,:)
      n1=n+1
      pi8=8.d0*pi
      spi8=sqrt(pi8)
      allocate(vp(n1,n1))
      allocate(d(n1))
c --------------------- Trace of PV --------------------------------
      s=pi8*q2aver*p00
      do i=1,n
        s=s+p(i,i)*v_opt_e(i)
      enddo
c --------------------- Trace of LOG(1-VP) -------------------------
      s=s+log(pi8)-2.d0*logq_aver
      do i=1,n
        s=s+log(v_opt_e(i))
      enddo
      vp(1,1)=-p00
      vp(1,2:n1)=-p01
      vp(2:n1,1)=-p01
      do j=2,n1
        do i=2,n1
          vp(i,j)=-p(i-1,j-1)
        enddo
        vp(j,j)=vp(j,j)+1.d0/v_opt_e(j-1)
      enddo
      call eig_val_solver_r(n1,n1,vp,d)
      do i=1,n1
        s=s+log(d(i))
      enddo
      deallocate(d,vp)
      s=s/2.d0
      phi_functional_c=phi_functional_c+et_sum_nu(0)*s*wgt(1)*2.d0
      end
