      subroutine boson_pack_tau(alf,f_cmpl,n0,bet,f_real,n1,n)
c     ------- Packing from P(t) into A(tau)  -------------------
c	  Up  triangle of A =   2 * real parts of P
c	  Low triangle of A =   2 * imaginary parts of P
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: alf,bet
	  complex*16, intent(in) :: f_cmpl(n0,n0)
	  real*8, intent(inout) :: f_real(n1,n1)
      integer :: i,j
      real*8 :: a
      a=2.d0*alf
      if(bet==0.d0) f_real=0.d0
	  do j=1,n
	    do i=1,j-1
	      f_real(i,j)=bet*f_real(i,j)+a*real(f_cmpl(i,j))
	      f_real(j,i)=bet*f_real(j,i)+a*imag(f_cmpl(i,j))
	    enddo
	    f_real(j,j)=bet*f_real(j,j)+a*real(f_cmpl(j,j))
	  enddo
      end




      subroutine boson_pack_tau3(alf,cmm,cmi,cii,n0,bet,rmm,rmi,rii,n1,
     &                           n)
c   ------- Packing from P(t) into A(tau)  -------------------
c	Up  triangle of A =   2 * real parts of P
c	Low triangle of A =   2 * imaginary parts of P
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: alf,bet
	complex*16, intent(in) :: cmm(n_pbmt,n_pbmt),cmi(n_pbmt,n0),
     &                          cii(n0,n0)
	real*8, intent(inout) :: rii(n1,n1)
	complex*16, intent(inout) :: rmm(ndim_pbmt(me_b+1)),
     &                             rmi(n_pbmt,n1)
      integer :: i,j,ij,ij0
      real*8 :: a
      a=2.d0*alf
      if(bet==0.d0) then
        rmm=(0.d0,0.d0)
        rmi=(0.d0,0.d0)
        rii=0.d0
      endif
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          ij0=ij-n_mpi_pbmt(me_b+1)
          if(ij0<=0) cycle
          if(ij0>ndim_pbmt(me_b+1)) cycle
	    rmm(ij0)=bet*rmm(ij0)+a*cmm(i,j)
	  enddo
	enddo
	do j=1,n
	  do i=1,n_pbmt
	    rmi(i,j)=bet*rmi(i,j)+a*cmi(i,j)
	  enddo
	enddo
	do j=1,n
	  do i=1,j-1
	    rii(i,j)=bet*rii(i,j)+a*real(cii(i,j))
	    rii(j,i)=bet*rii(j,i)+a*imag(cii(i,j))
	  enddo
	  rii(j,j)=bet*rii(j,j)+a*real(cii(j,j))
	enddo
      end
