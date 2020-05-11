      subroutine pack_hermit(a_compl,a_real,n,n0,n1,alpha,beta)
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: alpha,beta
	  complex*16, intent(in) :: a_compl(n0,n0)
	  real*8, intent(out) :: a_real(n1,n1)
	  integer :: i,j
	  if(alpha==0.d0) a_real=0.d0
	  do j=1,n
	    do i=1,j-1
	      a_real(i,j)=alpha*a_real(i,j)+beta*real(a_compl(i,j))
	      a_real(j,i)=alpha*a_real(j,i)+beta*imag(a_compl(i,j))
	    enddo
	    a_real(j,j)=alpha*a_real(j,j)+beta*real(a_compl(j,j))
	  enddo
      end
      
      
      
	  subroutine unpack_hermit(a_compl,a_real,n,n0,n1,alpha,beta)
      implicit none
	  integer, intent(in) :: n,n0,n1
	  complex*16, intent(in) :: alpha,beta
	  real*8, intent(in) :: a_real(n0,n0)
	  complex*16, intent(out) :: a_compl(n1,n1)
	  integer :: i,j
	  if(alpha==(0.d0,0.d0)) a_compl=(0.d0,0.d0)
	  do j=1,n
	    do i=1,j-1
	      a_compl(i,j)=alpha*a_compl(i,j)
     &		          +beta*dcmplx(a_real(i,j),a_real(j,i))
	      a_compl(j,i)=conjg(a_compl(i,j))
	    enddo
	    a_compl(j,j)=alpha*a_compl(j,j)+beta*a_real(j,j)
	  enddo
      end

      subroutine pack_sym(a_compl,a_real,n,n0,n1,alpha,beta)
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: alpha,beta
	real*8, intent(in) :: a_compl(n0,n0)
	real*8, intent(out) :: a_real(n1,n1)
	integer :: i,j
	if(alpha==0.d0) a_real=0.d0
	do j=1,n
	  do i=1,j
	    a_real(i,j)=alpha*a_real(i,j)+beta*a_compl(i,j)
	  enddo
	enddo
      end
      
      
      
	  subroutine unpack_hermit3(work,rmm,rmi,rii,n)
      use atom_mod
      use parallel_mod
      implicit none
	  integer, intent(in) :: n
      complex*16, intent(in) :: rmm(ndim_pbmt(me_b+1)),
     &                          rmi(nd_b_pbmt(me_b+1),n),
     &                          rii(nd_b_pbint(me_b+1))
	  complex*16, intent(out) :: work(n_pbmt+n,n_pbmt+n)
	  integer :: i,j,ij,ind_i,n1,n2
	  complex*16 :: cc
      work=(0.d0,0.d0)
c ------------------------------------------------------------------
      n1=n_mpi_pbmt(me_b+1)
      n2=n1+ndim_pbmt(me_b+1)
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          if(ij>n1.and.ij<=n2) then
            work(i,j)=rmm(ij-n1)
            if(j/=i) work(j,i)=conjg(rmm(ij-n1))
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      do j=1,n
        do ind_i=1,nd_b_pbmt(me_b+1)
          i=n0_b_pbmt(me_b+1)+ind_i
          work(i,n_pbmt+j)=rmi(ind_i,j)
          work(n_pbmt+j,i)=conjg(rmi(ind_i,j))
        enddo
      enddo
c ------------------------------------------------------------------
      n1=n0_b_pbint(me_b+1)
      n2=n1+nd_b_pbint(me_b+1)
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij>n1.and.ij<=n2) then
            cc=rii(ij-n1)
            work(n_pbmt+i,n_pbmt+j)=cc
            work(n_pbmt+j,n_pbmt+i)=conjg(cc)
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      if(nproc_b/=1) call dgop(work,2*(n_pbmt+n)**2,'  +',comm_b)
      end
      
      
      
      


      subroutine unpack_hermit3r(work,rmm,rmi,rii,n)
      use atom_mod
      use parallel_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: rmm(ndim_pbmt(me_b+1)),
     &                      rmi(nd_b_pbmt(me_b+1),n),
     &                      rii(nd_b_pbint(me_b+1))
      real*8, intent(out) :: work(n_pbmt+n,n_pbmt+n)
      integer :: i,j,ij,n1,n2,ind_i
      work=0.d0
c ------------------------------------------------------------------
      n1=n_mpi_pbmt(me_b+1)
      n2=n1+ndim_pbmt(me_b+1)
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          if(ij>n1.and.ij<=n2) then
            work(i,j)=rmm(ij-n1)
            if(j/=i) work(j,i)=rmm(ij-n1)
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      do j=1,n
        do ind_i=1,nd_b_pbmt(me_b+1)
          i=n0_b_pbmt(me_b+1)+ind_i
          work(i,n_pbmt+j)=rmi(ind_i,j)
          work(n_pbmt+j,i)=rmi(ind_i,j)
        enddo
      enddo
c ------------------------------------------------------------------
      n1=n0_b_pbint(me_b+1)
      n2=n1+nd_b_pbint(me_b+1)
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij>n1.and.ij<=n2) then
            work(n_pbmt+i,n_pbmt+j)=rii(ij-n1)
            work(n_pbmt+j,n_pbmt+i)=rii(ij-n1)
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      if(nproc_b/=1) call dgop(work,(n_pbmt+n)**2,'  +',comm_b)
      end

      subroutine pack_hermit3(a_compl,rmm,rmi,rii,n)
      use atom_mod
      use parallel_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: a_compl(n_pbmt+n,n_pbmt+n)
	  complex*16, intent(out) :: rmm(ndim_pbmt(me_b+1)),
     &                           rmi(nd_b_pbmt(me_b+1),n),
     &                           rii(nd_b_pbint(me_b+1))
	  integer :: i,j,i1,j1,ij,ij0,ind_i
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          ij0=ij-n_mpi_pbmt(me_b+1)
          if(ij0<=0) cycle
          if(ij0>ndim_pbmt(me_b+1)) cycle
	      rmm(ij0)=a_compl(i,j)
	    enddo
	  enddo
	  do j=1,n
        j1=n_pbmt+j
        do ind_i=1,nd_b_pbmt(me_b+1)
          i=n0_b_pbmt(me_b+1)+ind_i
	      rmi(ind_i,j)=a_compl(i,j1)
	    enddo
	  enddo
	  ij=0
	  do j=1,n
        j1=n_pbmt+j
	    do i=1,j
          i1=n_pbmt+i
          ij=ij+1
          ij0=ij-n0_b_pbint(me_b+1)
          if(ij0<=0) cycle
          if(ij0>nd_b_pbint(me_b+1)) cycle
	      rii(ij0)=a_compl(i1,j1)
	    enddo
	  enddo
      end



      subroutine pack_hermit3r(a_compl,rmm,rmi,rii,n)
      use atom_mod
      use parallel_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: a_compl(n_pbmt+n,n_pbmt+n)
      real*8, intent(out) :: rmm(ndim_pbmt(me_b+1)),
     &                       rmi(nd_b_pbmt(me_b+1),n),
     &                       rii(nd_b_pbint(me_b+1))
      integer :: i,j,i1,j1,ij,ij0,ind_i
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          ij0=ij-n_mpi_pbmt(me_b+1)
          if(ij0<=0) cycle
          if(ij0>ndim_pbmt(me_b+1)) cycle
          rmm(ij0)=a_compl(i,j)
        enddo
      enddo
      do j=1,n
        j1=n_pbmt+j
        do ind_i=1,nd_b_pbmt(me_b+1)
          i=n0_b_pbmt(me_b+1)+ind_i
          rmi(ind_i,j)=a_compl(i,j1)
        enddo
      enddo
      ij=0
      do j=1,n
        j1=n_pbmt+j
        do i=1,j
          i1=n_pbmt+i
          ij=ij+1
          ij0=ij-n0_b_pbint(me_b+1)
          if(ij0<=0) cycle
          if(ij0>nd_b_pbint(me_b+1)) cycle
          rii(ij0)=a_compl(i1,j1)
        enddo
      enddo
      end

      subroutine unpack_hermit3rOLD(a_compl,rmm,rmi,rii,n,iq)
      use atom_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: n,iq
      real*8, intent(in) :: rii(nd_b_pbint(me_b+1)),
     &                      rmi(nd_b_pbmt(me_b+1),n),
     &                      rmm(ndim_pbmt(me_b+1))
      complex*16, intent(out) :: a_compl(n_pbmt+n,n_pbmt+n)
      integer :: i,j,i1,j1,ind_ij,ij,n1,n2
      real*8, allocatable :: tr(:,:),tt(:)
      complex*16, allocatable :: t2(:,:)
      allocate(tt(n_pbmt*(n_pbmt+1)/2))
      tt=0.d0
      do ind_ij=1,ndim_pbmt(me_b+1)
        ij=n_mpi_pbmt(me_b+1)+ind_ij
        tt(ij)=rmm(ind_ij)
      enddo
      if(nproc_b/=1) call dgop(tt,n_pbmt*(n_pbmt+1)/2,'  +',comm_b)
      allocate(tr(n_pbmt,n_pbmt))
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          tr(i,j)=tt(ij)
          if(j/=i) tr(j,i)=tt(ij)
        enddo
      enddo
      deallocate(tt)
      allocate(t2(n_pbmt,n_pbmt))
      call pb_c_from_r_mm(t2,tr,pnt(1,iq))
      a_compl(1:n_pbmt,1:n_pbmt)=t2
      deallocate(tr,t2)
c ------------------------------------------------------------------
      allocate(t2(n_pbmt,n))
      call p_c_from_r_mi(n,t2,rmi,iq)
      do j=1,n
        j1=n_pbmt+j
        do i=1,n_pbmt
          a_compl(i,j1)=t2(i,j)
          a_compl(j1,i)=conjg(a_compl(i,j1))
        enddo
      enddo
      deallocate(t2)
c ------------------------------------------------------------------
      allocate(tr(n,n))
      n1=n0_b_pbint(me_b+1)
      n2=n1+nd_b_pbint(me_b+1)
      tr=0.d0
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij>n1.and.ij<=n2) then
            tr(i,j)=rii(ij-n1)
            tr(j,i)=rii(ij-n1)
          endif
        enddo
      enddo
      if(nproc_b/=1) call dgop(tr,n*n,'  +',comm_b)
      ij=0
      do j=1,n
        j1=n_pbmt+j
        do i=1,n
          i1=n_pbmt+i
          a_compl(i1,j1)=tr(i,j)
        enddo
      enddo
      deallocate(tr)
      end
