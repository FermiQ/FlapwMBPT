      subroutine phi_rpa_c_sclp(n,nrow,ncol,k,i_nu,p,v,desc,refrow,
     &                          refcol)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol,i_nu,k,desc(9),refrow(nrow),
     &                       refcol(ncol)
      complex*16, intent(in) :: p(nrow,ncol),v(nrow,ncol)
      integer :: i,j,ii,jj
      real*8 :: s
      real*8, allocatable :: d(:)
      complex*16, allocatable :: pc(:,:),vp(:,:),a(:,:)
      allocate(pc(nrow,ncol))
      allocate(vp(nrow,ncol))
      allocate(a(nrow,ncol))
      allocate(d(n))
#ifdef SCLP
      call pzgemm('n','n',n,n,n,(1.d0,0.d0),v,1,1,desc,p,1,1,desc,
     &            (0.d0,0.d0),vp,1,1,desc)
      call pzgemm('n','n',n,n,n,(1.d0,0.d0),p,1,1,desc,v,1,1,desc,
     &            (0.d0,0.d0),pc,1,1,desc)
#endif
      s=0.d0
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          if(ii==jj) s=s+vp(i,j)
        enddo
      enddo
      call dgop(s,1,'  +',comm_b)
#ifdef SCLP
      call pzgemm('n','n',n,n,n,(1.d0,0.d0),pc,1,1,desc,vp,1,1,desc,
     &            (0.d0,0.d0),a,1,1,desc)
#endif
      a=a-vp-pc
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          if(ii==jj) a(i,j)=(1.d0,0.d0)+a(i,j)
        enddo
      enddo
#ifdef SCLP
      call eig_val_only_sclp(n,a,d,nrow,ncol,desc)
#endif
      d=log(d)
      s=s/2.d0+sum(d)/4.d0
      phi_functional_c=phi_functional_c
     &                +et_sum_nu(i_nu)*s*wgt(k)*2.d0
      deallocate(pc,vp,a,d)
      end


      subroutine phi_rpa_c_sclp_r(n,nrow,ncol,k,i_nu,p,v,desc,refrow,
     &                            refcol)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol,i_nu,k,desc(9),refrow(nrow),
     &                       refcol(ncol)
      real*8, intent(in) :: p(nrow,ncol),v(nrow,ncol)
      integer :: i,j,ii,jj
      real*8 :: s
      real*8, allocatable :: d(:),pc(:,:),vp(:,:),a(:,:)
      allocate(pc(nrow,ncol))
      allocate(vp(nrow,ncol))
      allocate(a(nrow,ncol))
      allocate(d(n))
#ifdef SCLP
      call pdgemm('n','n',n,n,n,1.d0,v,1,1,desc,p,1,1,desc,0.d0,vp,1,1,
     &            desc)
      call pdgemm('n','n',n,n,n,1.d0,p,1,1,desc,v,1,1,desc,0.d0,pc,1,1,
     &            desc)
#endif
      s=0.d0
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          if(ii==jj) s=s+vp(i,j)
        enddo
      enddo
      call dgop(s,1,'  +',comm_b)
#ifdef SCLP
      call pdgemm('n','n',n,n,n,1.d0,pc,1,1,desc,vp,1,1,desc,0.d0,a,1,1,
     &            desc)
#endif
      a=a-vp-pc
      do j=1,ncol
        jj=refcol(j)
        do i=1,nrow
          ii=refrow(i)
          if(ii==jj) a(i,j)=1.d0+a(i,j)
        enddo
      enddo
#ifdef SCLP
      call eig_val_only_sclp_r(n,a,d,nrow,ncol,desc)
#endif
      d=log(d)
      s=s/2.d0+sum(d)/4.d0
      phi_functional_c=phi_functional_c
     &                +et_sum_nu(i_nu)*s*wgt(k)*2.d0
      deallocate(pc,vp,a,d)
      end
