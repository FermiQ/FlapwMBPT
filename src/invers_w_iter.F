      subroutine invers_w_iter(n,a,eps,eps1,itmax,nit,nn,n0,nd)
c     Replaces matrix A with (1-A)^(-1)
c     Matrix A is real symmetric
c     Parameter EPS defines the importance of matrix elements through
c     the tolerance: TOL=EPS*maxval(ABS(A))
c     EPS1 defines the convergence
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,itmax,n0,nd
      real*8, intent(in) :: eps,eps1
      real*8, intent(inout) :: a(nd,n)
      integer, intent(out) :: nit,nn
      integer :: i,j,it,i0,i1,j0,j1,n1,ind_i,ib,ie,nni
      real*8 :: tol,t,da
      integer, allocatable :: flag(:),perm(:),perm1(:)
      real*8,allocatable :: at(:,:),dd(:,:),
     &                      ax21(:,:),ax2(:,:),adiag(:),
     &                      t1(:,:),t2(:,:),d(:)
      ib=n0+1
      ie=n0+nd
      tol=eps*maxval(abs(a))
      if(nproc_b/=1) call dgop(tol,1,'max',comm_b)
      allocate(flag(n))
      flag=0
      do j=1,n
        do ind_i=1,nd
          i=n0+ind_i
          if(i/=j.and.abs(a(ind_i,j))>tol) then
            flag(i)=flag(i)+1
            flag(j)=flag(j)+1
          endif
        enddo
      enddo
      if(nproc_b/=1) call igop(flag,n,'  +',comm_b)
      nn=0
      do i=1,n
        if(flag(i)/=0) nn=nn+1
      enddo
      allocate(perm(n))
      j=0
      do i=1,n
        if(flag(i)/=0) then
          j=j+1
          perm(j)=i
        endif
      enddo
      do i=1,n
        if(flag(i)==0) then
          j=j+1
          perm(j)=i
        endif
      enddo
      deallocate(flag)
c ------------------------------------------------------------------
      n1=n-nn
      allocate(at(nn,nn))
      at=0.d0
      do j=1,nn
        j0=perm(j)
        do i=1,nn
          i0=perm(i)
          if(i0>=ib.and.i0<=ie) at(i,j)=a(i0-n0,j0)
        enddo
      enddo
      if(nproc_b/=1) call dgop(at,nn*nn,'  +',comm_b)
      at=-at
      do i=1,nn
        at(i,i)=at(i,i)+1.d0
      enddo
      call invers_r(nn,at,nn)
      allocate(adiag(n1))
      adiag=0.d0
      do i=1,n1
        i0=perm(nn+i)
        if(i0>=ib.and.i0<=ie) adiag(i)=1.d0/(1.d0-a(i0-n0,i0))
      enddo
      if(nproc_b/=1) call dgop(adiag,n1,'  +',comm_b)
c ------------------------------------------------------------------
      allocate(ax21(n1,nn))
      allocate(ax2(n,n1))
      ax21=0.d0
      do j=1,nn
        j1=perm(j)
        do i=1,n1
          i1=perm(nn+i)
          if(i1>=ib.and.i1<=ie) ax21(i,j)=a(i1-n0,j1)
        enddo
      enddo
      if(nproc_b/=1) call dgop(ax21,n1*nn,'  +',comm_b)
      ax2=0.d0
      do j=1,n1
        j1=perm(nn+j)
        do i=1,nn
          i1=perm(i)
          if(i1>=ib.and.i1<=ie) ax2(i,j)=a(i1-n0,j1)
        enddo
      enddo
      do j=1,n1
        j1=perm(nn+j)
        do i=nn+1,n
          if(i==nn+j) then
            ax2(i,j)=0.d0
          else
            i1=perm(i)
            if(i1>=ib.and.i1<=ie) ax2(i,j)=a(i1-n0,j1)
          endif
        enddo
      enddo
      if(nproc_b/=1) call dgop(ax2,n1*n,'  +',comm_b)
c ------------------------------------------------------------------
      allocate(t1(n1,nn))
      call dgemm('n','n',n1,nn,nn,1.d0,ax21,n1,at,nn,0.d0,t1,n1)
      ax21=t1
      deallocate(t1)
      do j=1,n1
        t=adiag(j)
        ax2(:,j)=ax2(:,j)*t
      enddo
c ---- First (fast) iteration --------------------------------------
      allocate(dd(nd,n))
      dd=0.d0
      do i=1,n1
        nni=nn+i
        if(nni>=ib.and.nni<=ie) then
          t=adiag(i)
          dd(nni-n0,1:nn)=t*ax21(i,:)
          dd(nni-n0,nn+1:n)=t*ax2(nni,:)
        endif
      enddo
      if(ib<=nn) then
        i1=min(nd,nn-n0)
        call dgemm('n','n',i1,n1,nn,1.d0,at(ib,1),nn,ax2,n,0.d0,
     &             dd(1,nn+1),nd)
        a(1:i1,1:nn)=at(ib:ib+i1-1,:)
        a(1:i1,nn+1:n)=dd(1:i1,nn+1:n)
      endif
      if(ie>nn) then
        i1=1
        if(ib<=nn) i1=nn+1-n0
        a(i1:nd,1:n)=dd(i1:nd,1:n)
      endif
      do i=1,n1
        nni=nn+i
        if(nni>=ib.and.nni<=ie) a(nni-n0,nni)=a(nni-n0,nni)+adiag(i)
      enddo
      da=maxval(abs(dd))
      if(nproc_b/=1) call dgop(da,1,'max',comm_b)
      nit=1
      deallocate(adiag)
      if(da<eps1) goto 1
c ----- The rest of iterations -------------------------------------
      allocate(t1(nd,nn))
      allocate(t2(nd,n1))
      do it=2,itmax
        nit=nit+1
        call dgemm('n','n',nd,nn,n1,1.d0,dd(1,nn+1),nd,ax21,n1,0.d0,t1,
     &             nd)
        call dgemm('n','n',nd,n1,n,1.d0,dd,nd,ax2,n,0.d0,t2,nd)
        dd(:,1:nn)=t1
        dd(:,nn+1:n)=t2
        a=a+dd
        da=maxval(abs(dd))
        if(nproc_b/=1) call dgop(da,1,'max',comm_b)
        if(da<eps1) goto 1
      enddo
      deallocate(t1,t2)
c ------------------------------------------------------------------
1     continue
      deallocate(at,ax21,ax2,dd)
c ----- Backwards permutation of the matrix A ----------------------
      allocate(perm1(n))
      do i=1,n
        j=perm(i)
        perm1(j)=i
      enddo
      deallocate(perm)
      allocate(d(n))
      do j=1,n
        d=0.d0
        d(ib:ie)=a(:,j)
        if(nproc_b/=1) call dgop(d,n,'  +',comm_b)
        do ind_i=1,nd
          i=n0+ind_i
          i1=perm1(i)
          a(ind_i,j)=d(i1)
        enddo
      enddo
      do j=1,nd
        d=a(j,:)
        do i=1,n
          i1=perm1(i)
          a(j,i)=d(i1)
        enddo
      enddo
      deallocate(perm1,d)
      end
