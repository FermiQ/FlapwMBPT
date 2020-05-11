      subroutine boson_unpack_tau(f_cmpl,n1,f_real,n0,n)
c   ------- Unpacking from A(t) into P(t)  -------------------
c	  Up  triangle of P_real =    2 * real parts of P
c	  Low triangle of P_real =    2 * imaginary parts of P
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f_real(n0,n0)
	  complex*16, intent(out) :: f_cmpl(n1,n1)
      integer :: i,j
	  do j=1,n
	    do i=1,j-1
	      f_cmpl(i,j)=0.5d0*dcmplx(f_real(i,j),f_real(j,i))
	      f_cmpl(j,i)=conjg(f_cmpl(i,j))
	    enddo
	    f_cmpl(j,j)=0.5d0*f_real(j,j)
	  enddo
      end




      subroutine boson_unpack_tau3(cmm,cmi,cii,n1,rmm,rmi,rii,n0,n)
c   ------- Unpacking from A(t) into P(t)  -------------------
c	  Up  triangle of P_real =    2 * real parts of P
c	  Low triangle of P_real =    2 * imaginary parts of P
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
      complex*16, intent(in) :: rmi(nd_b_pbmt(me_b+1),n0),
     &                          rii(nd_b_pbint(me_b+1)),
     &                          rmm(ndim_pbmt(me_b+1))
	  complex*16, intent(out) :: cmm(n_pbmt,n_pbmt),cmi(n_pbmt,n1),
     &                           cii(n1,n1)
      integer :: i,j,ij,ind_i,n2,n11
      complex*16 :: cc
      n11=n_mpi_pbmt(me_b+1)
      n2=n11+ndim_pbmt(me_b+1)
      cmm=(0.d0,0.d0)
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          if(ij>n11.and.ij<=n2) then
            cmm(i,j)=rmm(ij-n11)/2
            if(j/=i) cmm(j,i)=conjg(rmm(ij-n11))/2
          endif
        enddo
      enddo
      if(nproc_b/=1) call dgop(cmm,2*n_pbmt**2,'  +',comm_b)
c ------------------------------------------------------------------
      cmi=(0.d0,0.d0)
      do j=1,n
        do ind_i=1,nd_b_pbmt(me_b+1)
          i=n0_b_pbmt(me_b+1)+ind_i
          cmi(i,j)=rmi(ind_i,j)/2
        enddo
      enddo
      if(nproc_b/=1) call dgop(cmi,2*n_pbmt*n1,'  +',comm_b)
c ------------------------------------------------------------------
      n11=n0_b_pbint(me_b+1)
      n2=n11+nd_b_pbint(me_b+1)
      cii=(0.d0,0.d0)
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij>n11.and.ij<=n2) then
            cc=0.5d0*rii(ij-n11)
            cii(i,j)=cc
            cii(j,i)=conjg(cc)
          endif
        enddo
      enddo
      if(nproc_b/=1) call dgop(cii,2*n1*n1,'  +',comm_b)
      end




      subroutine boson_unpack_tau3all(work,rmm,rmi,rii,n)
c   ------- Unpacking from A(t) into P(t)  -------------------
c     Up  triangle of P_real =    2 * real parts of P
c     Low triangle of P_real =    2 * imaginary parts of P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      complex*16, intent(in) :: rmi(nd_b_pbmt(me_b+1),nplw_gw),
     &                          rmm(ndim_pbmt(me_b+1)),
     &                          rii(nd_b_pbint(me_b+1))
      complex*16, intent(out) :: work(n_pbmt+n,n_pbmt+n)
      integer :: i,j,ij,n1,n2,ind_i
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
            work(i,j)=rmm(ij-n1)/2
            if(j/=i) work(j,i)=conjg(rmm(ij-n1))/2
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      do j=1,n
        do ind_i=1,nd_b_pbmt(me_b+1)
          i=n0_b_pbmt(me_b+1)+ind_i
          work(i,n_pbmt+j)=rmi(ind_i,j)/2
          work(n_pbmt+j,i)=conjg(rmi(ind_i,j))/2
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
            cc=0.5d0*rii(ij-n1)
            work(n_pbmt+i,n_pbmt+j)=cc
            work(n_pbmt+j,n_pbmt+i)=conjg(cc)
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      if(nproc_b/=1) call dgop(work,2*(n_pbmt+n)**2,'  +',comm_b)
      end






      subroutine boson_unpack_tau3allr(work,rmm,rmi,rii,n)
c   ------- Unpacking from A(t) into P(t)  -------------------
c     Up  triangle of P_real =    2 * real parts of P
c     Low triangle of P_real =    2 * imaginary parts of P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: rmi(nd_b_pbmt(me_b+1),nplw_gw),
     &                      rmm(ndim_pbmt(me_b+1)),
     &                      rii(nd_b_pbint(me_b+1))
      real*8, intent(out) :: work(n_pbmt+n,n_pbmt+n)
      integer :: i,j,ind_i,ij,n1,n2
      real*8 :: cc
      work=0.d0
c ------------------------------------------------------------------
      n1=n_mpi_pbmt(me_b+1)
      n2=n1+ndim_pbmt(me_b+1)
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          if(ij>n1.and.ij<=n2) then
            work(i,j)=rmm(ij-n1)/2
            if(j/=i) work(j,i)=rmm(ij-n1)/2
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      do j=1,n
        do ind_i=1,nd_b_pbmt(me_b+1)
          i=n0_b_pbmt(me_b+1)+ind_i
          work(i,n_pbmt+j)=rmi(ind_i,j)/2
          work(n_pbmt+j,i)=rmi(ind_i,j)/2
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
            cc=0.5d0*rii(ij-n1)
            work(n_pbmt+i,n_pbmt+j)=cc
            work(n_pbmt+j,n_pbmt+i)=cc
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      if(nproc_b/=1) call dgop(work,(n_pbmt+n)**2,'  +',comm_b)
      end




      subroutine boson_unpack_tau3r(cmm,cmi,cii,n1,rmm,rmi,rii,n0,n,iq)
c   ------- Unpacking from A(t) into P(t)  -------------------
c     Up  triangle of P_real =    2 * real parts of P
c     Low triangle of P_real =    2 * imaginary parts of P
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0,n1,iq
      real*8, intent(in) :: rmi(nd_b_pbmt(me_b+1),n0),
     &                      rmm(ndim_pbmt(me_b+1)),
     &                      rii(nd_b_pbint(me_b+1))
      complex*16, intent(out) :: cmm(n_pbmt,n_pbmt),cmi(n_pbmt,n1),
     &                           cii(n1,n1)
      integer :: i,j,ij,n11,n2
      complex*16 :: cc
      real*8, allocatable :: tr(:,:)
      n11=n_mpi_pbmt(me_b+1)
      n2=n11+ndim_pbmt(me_b+1)
      allocate(tr(n_pbmt,n_pbmt))
      tr=0.d0
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          if(ij>n11.and.ij<=n2) then
            tr(i,j)=rmm(ij-n11)/2
            if(j/=i) tr(j,i)=rmm(ij-n11)/2
          endif
        enddo
      enddo
      if(nproc_b/=1) call dgop(tr,n_pbmt**2,'  +',comm_b)
      call pb_c_from_r_mm(cmm,tr,pnt(1,iq))
      deallocate(tr)
c ------------------------------------------------------------------
      call p_c_from_r_mi(n,cmi,rmi,iq)
      cmi=cmi/2
c ------------------------------------------------------------------
      n11=n0_b_pbint(me_b+1)
      n2=n11+nd_b_pbint(me_b+1)
      cii=0.d0
      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          if(ij>n11.and.ij<=n2) then
            cc=0.5d0*rii(ij-n11)
            cii(i,j)=cc
            cii(j,i)=cc
          endif
        enddo
      enddo
      if(nproc_b/=1) call dgop(cii,2*n1*n1,'  +',comm_b)
      end
