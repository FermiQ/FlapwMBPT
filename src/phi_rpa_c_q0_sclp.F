      subroutine phi_rpa_c_q0_sclp(n,i_nu,p,p_head,p_wing,nrow,ncol,
     &                             refrow,refcol)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,i_nu,nrow,ncol,refrow(nrow),refcol(ncol)
      complex*16, intent(in) :: p(nrow,ncol),p_head(3,3),p_wing(3,n)
      integer :: i,j,n1,ii,jj
      real*8 :: s,p00,pi8,spi8
      real*8, allocatable :: d(:),sv(:)
      complex*16, allocatable :: vp(:,:)
      n1=n+1
      pi8=8.d0*pi
      spi8=sqrt(pi8)
      p00=0.d0
      do i=1,3
        p00=p00+p_head(i,i)
      enddo
      p00=p00/3.d0
      allocate(vp(n1,n1))
      allocate(d(n1))
      allocate(sv(n))
      do i=1,n
        sv(i)=sqrt(v_opt_e(i))
      enddo
      vp=(0.d0,0.d0)
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          vp(ii+1,jj+1)=sv(ii)*p(i,j)*sv(jj)
        enddo
      enddo
      call dgop(vp,2*n1*n1,'  +',comm_b)
      vp(1,1)=pi8*p00
      do i=2,n1
        vp(1,i)=spi8*p_wing(3,i-1)*sv(i-1)
        vp(i,1)=conjg(vp(1,i))
      enddo
c --------------------- Trace of PV --------------------------------
      s=0.d0
      do i=1,n1
        s=s+vp(i,i)
      enddo
c --------------------- Trace of LOG(1-VP) -------------------------
      vp=-vp
      do i=1,n1
        vp(i,i)=(1.d0,0.d0)+vp(i,i)
      enddo
      call eig_val_solver(n1,n1,vp,d)
      do i=1,n1
        s=s+log(d(i))
      enddo
      deallocate(sv)
      deallocate(d,vp)
      s=s/2.d0
      phi_functional_c=phi_functional_c
     &                +et_sum_nu(i_nu)*s*wgt(1)*2.d0
      end


      subroutine phi_rpa_c_q0_sclp_r(n,i_nu,p,p_head,p_wing,nrow,ncol,
     &                               refrow,refcol)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,i_nu,nrow,ncol,refrow(nrow),refcol(ncol)
      real*8, intent(in) :: p(nrow,ncol),p_head(3,3),p_wing(3,n)
      integer :: i,j,n1,ii,jj
      real*8 :: s,p00,pi8,spi8
      real*8, allocatable :: d(:),sv(:),vp(:,:)
      n1=n+1
      pi8=8.d0*pi
      spi8=sqrt(pi8)
      p00=0.d0
      do i=1,3
        p00=p00+p_head(i,i)
      enddo
      p00=p00/3.d0
      allocate(vp(n1,n1))
      allocate(d(n1))
      allocate(sv(n))
      do i=1,n
        sv(i)=sqrt(v_opt_e(i))
      enddo
      vp=0.d0
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          vp(ii+1,jj+1)=sv(ii)*p(i,j)*sv(jj)
        enddo
      enddo
      call dgop(vp,n1*n1,'  +',comm_b)
      vp(1,1)=pi8*p00
      do i=2,n1
        vp(1,i)=spi8*p_wing(3,i-1)*sv(i-1)
        vp(i,1)=vp(1,i)
      enddo
c --------------------- Trace of PV --------------------------------
      s=0.d0
      do i=1,n1
        s=s+vp(i,i)
      enddo
c --------------------- Trace of LOG(1-VP) -------------------------
      vp=-vp
      do i=1,n1
        vp(i,i)=1.d0+vp(i,i)
      enddo
      call eig_val_solver_r(n1,n1,vp,d)
      do i=1,n1
        s=s+log(d(i))
      enddo
      deallocate(sv)
      deallocate(d,vp)
      s=s/2.d0
      phi_functional_c=phi_functional_c
     &                +et_sum_nu(i_nu)*s*wgt(1)*2.d0
      end
