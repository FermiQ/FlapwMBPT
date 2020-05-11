      subroutine phi_rpa_c_q0_met_sclp(n,p,p00,p01,nrow,ncol,refrow,
     &                                 refcol)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol,refrow(nrow),refcol(ncol)
      real*8, intent(in) :: p00
      complex*16, intent(in) :: p(nrow,ncol),p01(n)
      integer :: i,j,n1,ii,jj
      real*8 :: s,pi8,spi8
      real*8, allocatable :: d(:)
      complex*16, allocatable :: vp(:,:)
      n1=n+1
      pi8=8.d0*pi
      spi8=sqrt(pi8)
      allocate(vp(n1,n1))
      allocate(d(n1))
c --------------------- Trace of PV --------------------------------
      s=0.d0
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          if(ii==jj) s=s+p(i,j)*v_opt_e(ii)
        enddo
      enddo
      call dgop(s,1,'  +',comm_b)
      s=s+pi8*q2aver*p00
c --------------------- Trace of LOG(1-VP) -------------------------
      s=s+log(pi8)-2.d0*logq_aver
      do i=1,n
        s=s+log(v_opt_e(i))
      enddo
      vp=(0.d0,0.d0)
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          vp(ii+1,jj+1)=-p(i,j)
        enddo
      enddo
      call dgop(vp,2*n1*n1,'  +',comm_b)
      do j=2,n1
        vp(j,j)=vp(j,j)+1.d0/v_opt_e(j-1)
      enddo
      vp(1,1)=-p00
      vp(1,2:n1)=-p01
      vp(2:n1,1)=-conjg(p01)
      call eig_val_solver(n1,n1,vp,d)
      do i=1,n1
        s=s+log(d(i))
      enddo
      deallocate(d,vp)
      s=s/2.d0
      phi_functional_c=phi_functional_c+et_sum_nu(0)*s*wgt(1)*2.d0
      end


      subroutine phi_rpa_c_q0_met_sclp_r(n,p,p00,p01,nrow,ncol,refrow,
     &                                   refcol)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol,refrow(nrow),refcol(ncol)
      real*8, intent(in) :: p00,p(nrow,ncol),p01(n)
      integer :: i,j,n1,ii,jj
      real*8 :: s,pi8,spi8
      real*8, allocatable :: d(:),vp(:,:)
      n1=n+1
      pi8=8.d0*pi
      spi8=sqrt(pi8)
      allocate(vp(n1,n1))
      allocate(d(n1))
c --------------------- Trace of PV --------------------------------
      s=0.d0
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          if(ii==jj) s=s+p(i,j)*v_opt_e(ii)
        enddo
      enddo
      call dgop(s,1,'  +',comm_b)
      s=s+pi8*q2aver*p00
c --------------------- Trace of LOG(1-VP) -------------------------
      s=s+log(pi8)-2.d0*logq_aver
      do i=1,n
        s=s+log(v_opt_e(i))
      enddo
      vp=0.d0
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          vp(ii+1,jj+1)=-p(i,j)
        enddo
      enddo
      call dgop(vp,n1*n1,'  +',comm_b)
      do j=2,n1
        vp(j,j)=vp(j,j)+1.d0/v_opt_e(j-1)
      enddo
      vp(1,1)=-p00
      vp(1,2:n1)=-p01
      vp(2:n1,1)=-p01
      call eig_val_solver_r(n1,n1,vp,d)
      do i=1,n1
        s=s+log(d(i))
      enddo
      deallocate(d,vp)
      s=s/2.d0
      phi_functional_c=phi_functional_c+et_sum_nu(0)*s*wgt(1)*2.d0
      end
