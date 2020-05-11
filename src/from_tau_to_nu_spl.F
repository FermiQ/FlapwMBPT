      subroutine tau_to_nu_spl(p_nu,p_tau,dp_tau,n,key)
c     key = 0  - Free boundary conditions
c     key = 2  - First order derivatives at boundaries are given
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,key
	  real*8, intent(in) :: p_tau(n,0:n_tau/2),dp_tau(n)
	  real*8, intent(out) :: p_nu(n,0:n_nu)
      real*8, allocatable :: y(:,:)
      allocate(y(n,0:n_tau/2+key))
      y(:,0:n_tau/2)=2.d0*p_tau
      if(key==2) then
        y(:,n_tau/2+1)=2.d0*dp_tau
        y(:,n_tau/2+2)=0.d0
      endif
      call dgemm('n','n',n,n_nu+1,n_tau/2+1+key,1.d0,y,n,
     &           ab_nu_from_tau(0,0,1),n_tau/2+1+key,0.d0,p_nu,n)
      deallocate(y)
      end



      subroutine tau_to_nu_gen_spl(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation ----------------------
c   ---- We transform totally general complex F_TAU with N elements ---
c	  F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [0:n_nu] form
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_tau(n,2,ndim3_tau)
	  complex*16, intent(out) :: f_nu(n,ndim3_nu)
      integer :: i_tau,ind_tau,ind_nu,i,i_nu
	  real*8, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	  allocate(at(n,ndim3_tau))
	  allocate(bt(n,ndim3_tau))
	  allocate(an(n,0:n_nu))
	  allocate(bn(n,0:n_nu))
c -------- Real part transformation -----------------------------
      do ind_tau=1,ndim3_tau
        i_tau=me_t*ndim3_tau+ind_tau-1
        at(:,ind_tau)=real(f_tau(:,1,ind_tau))+real(f_tau(:,2,ind_tau))
        bt(:,ind_tau)=real(f_tau(:,1,ind_tau))-real(f_tau(:,2,ind_tau))
	  enddo
      call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,at,n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,1),
     &           n_tau/2+1+dp_tau_to_nu,0.d0,an,n)
      call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,bt,n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,2),
     &           n_tau/2+1+dp_tau_to_nu,0.d0,bn,n)
      if(nproc_t/=1) then
        call dgop(an,n*(n_nu+1),'  +',comm_t)
        call dgop(bn,n*(n_nu+1),'  +',comm_t)
      endif
	  do i=1,n
	    do ind_nu=1,ndim3_nu
	      i_nu=me_t*ndim3_nu+ind_nu-1
	      f_nu(i,ind_nu)=dcmplx(an(i,i_nu),bn(i,i_nu))
	    enddo
	  enddo
c -------- Imaginary part transformation -----------------------------
      do ind_tau=1,ndim3_tau
        i_tau=me_t*ndim3_tau+ind_tau-1
        at(:,ind_tau)=imag(f_tau(:,1,ind_tau))+imag(f_tau(:,2,ind_tau))
        bt(:,ind_tau)=imag(f_tau(:,1,ind_tau))-imag(f_tau(:,2,ind_tau))
	  enddo
      call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,at,n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,1),
     &           n_tau/2+1+dp_tau_to_nu,0.d0,an,n)
      call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,bt,n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,2),
     &           n_tau/2+1+dp_tau_to_nu,0.d0,bn,n)
      if(nproc_t/=1) then
        call dgop(an,n*(n_nu+1),'  +',comm_t)
        call dgop(bn,n*(n_nu+1),'  +',comm_t)
      endif
	  do i=1,n
	    do ind_nu=1,ndim3_nu
	      i_nu=me_t*ndim3_nu+ind_nu-1
	      f_nu(i,ind_nu)=f_nu(i,ind_nu)+(0.d0,1.d0)
     &	                             *dcmplx(an(i,i_nu),bn(i,i_nu))
	    enddo
	  enddo
	  deallocate(at,bt,an,bn)
      end



      subroutine tau_to_nu_ba(f_nu,n0,f_tau,n1,n)
c   ------- Transform TAU ---> NU representation -------------------
c	  F_tau - is stored in A form
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f_tau(n1,n1,ndim3_tau)
	  real*8, intent(out) :: f_nu(n0,n0,ndim3_nu)
      integer :: ind_tau,ind_nu,j,i_nu
	  real*8, allocatable :: at(:,:),an(:,:)
	  allocate(at(n,ndim3_tau))
	  allocate(an(n,0:n_nu))
      do j=1,n
        do ind_tau=1,ndim3_tau
	      at(:,ind_tau)=f_tau(1:n,j,ind_tau)
	    enddo
        call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,at,n,
     &             ab_nu_from_tau(me_t*ndim3_tau,0,1),
     &             n_tau/2+1+dp_tau_to_nu,0.d0,an,n)
        if(nproc_t/=1) call dgop(an,n*(n_nu+1),'  +',comm_t)
        do ind_nu=1,ndim3_nu
	      i_nu=me_t*ndim3_nu+ind_nu-1
	      f_nu(1:n,j,ind_nu)=an(:,i_nu)
	    enddo
      enddo
	  deallocate(at,an)
      end



      subroutine tau_to_nu_bb(f_nu,n0,f_tau,n1,n)
c   ------- Transform TAU ---> NU representation ----------------------
c	  F_tau - is stored in B form
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f_tau(n1,n1,ndim3_tau)
	  real*8, intent(out) :: f_nu(n0,n0,ndim3_nu)
      integer :: ind_tau,ind_nu,j,i_nu
	  real*8, allocatable :: at(:,:),an(:,:)
	  allocate(at(n,ndim3_tau))
	  allocate(an(n,0:n_nu))
      do j=1,n
        do ind_tau=1,ndim3_tau
	      at(:,ind_tau)=f_tau(1:n,j,ind_tau)
	    enddo
        call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,at,n,
     &             ab_nu_from_tau(me_t*ndim3_tau,0,2),
     &             n_tau/2+1+dp_tau_to_nu,0.d0,an,n)
        if(nproc_t/=1) call dgop(an,n*(n_nu+1),'  +',comm_t)
        do ind_nu=1,ndim3_nu
	      i_nu=me_t*ndim3_nu+ind_nu-1
	      f_nu(1:n,j,ind_nu)=an(:,i_nu)
	    enddo
      enddo
	  deallocate(at,an)
      end



      subroutine tau_to_nu_full(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation ----------------------
c   ---- We transform totally general complex F_TAU with N elements ---
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [-n_nu(1) : n_nu(2)] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_nu(n,ndim3_nu,2)
      integer :: i_tau,ind_tau,ind_nu,i,i_nu
	complex*16, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(bt(n,ndim3_tau))
	allocate(an(n,0:n_nu))
	allocate(bn(n,0:n_nu))
      do ind_tau=1,ndim3_tau
        i_tau=me_t*ndim3_tau+ind_tau-1
	  at(:,ind_tau)=f_tau(:,1,ind_tau)+f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	enddo
      call dgemm('n','n',2*n,n_nu+1,ndim3_tau,1.d0,at,2*n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,1),
     &           n_tau/2+1+dp_tau_to_nu,0.d0,an,2*n)
      call dgemm('n','n',2*n,n_nu+1,ndim3_tau,1.d0,bt,2*n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,2),
     &           n_tau/2+1+dp_tau_to_nu,0.d0,bn,2*n)
      if(nproc_t/=1) then
        call dgop(an,2*n*(n_nu+1),'  +',comm_t)
        call dgop(bn,2*n*(n_nu+1),'  +',comm_t)
      endif
	do i=1,n
	  do ind_nu=1,ndim3_nu
	    i_nu=me_t*ndim3_nu+ind_nu-1
	    f_nu(i,ind_nu,1)=an(i,i_nu)-(0.d0,1.d0)*bn(i,i_nu)
	    f_nu(i,ind_nu,2)=an(i,i_nu)+(0.d0,1.d0)*bn(i,i_nu)
	  enddo
	enddo
	deallocate(at,bt,an,bn)
      end
      
      subroutine tau_to_nu_full_seq(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation --------------------------
c   ---- We transform totally general complex F_TAU with N elements ------
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [-n_nu(1) : n_nu(2)] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_tau(n,0:n_tau)
	complex*16, intent(out) :: f_nu(n,0:n_nu,2)
      integer :: i_tau,i,i_nu,i1_tau
	complex*16, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	allocate(at(n,0:n_tau/2))
	allocate(bt(n,0:n_tau/2))
	allocate(an(n,0:n_nu))
	allocate(bn(n,0:n_nu))
      do i_tau=0,n_tau/2
        i1_tau=n_tau-i_tau
	  at(:,i_tau)=f_tau(:,i_tau)+f_tau(:,i1_tau)
	  bt(:,i_tau)=f_tau(:,i_tau)-f_tau(:,i1_tau)
	enddo
	call dgemm('n','n',2*n,n_nu+1,n_tau/2+1,1.d0,at,2*n,
     &           ab_nu_from_tau(0,0,1),n_tau/2+1,0.d0,an,2*n)
	call dgemm('n','n',2*n,n_nu+1,n_tau/2+1,1.d0,bt,2*n,
     &           ab_nu_from_tau(0,0,2),n_tau/2+1,0.d0,bn,2*n)
	f_nu=(0.d0,0.d0)
	do i=1,n
	  do i_nu=0,n_nu
	    f_nu(i,i_nu,1)=an(i,i_nu)-(0.d0,1.d0)*bn(i,i_nu)
	    f_nu(i,i_nu,2)=an(i,i_nu)+(0.d0,1.d0)*bn(i,i_nu)
	  enddo
	enddo
	deallocate(at,bt,an,bn)
      end



      subroutine tau_to_nu_full_kk4(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation ----------------------
c   ---- We transform totally general complex F_TAU with N elements ---
c	  F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [-n_nu(1) : n_nu(2)] form
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_tau(n,2,0:n_tau/2)
	  complex*16, intent(out) :: f_nu(n,0:n_nu,2)
      integer :: i_tau,i,i_nu
	  complex*16, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	  allocate(at(n,0:n_tau/2))
	  allocate(bt(n,0:n_tau/2))
	  allocate(an(n,0:n_nu))
	  allocate(bn(n,0:n_nu))
      an=(0.d0,0.d0)
      bn=(0.d0,0.d0)
      do i_tau=0,n_tau/2
	    at(:,i_tau)=f_tau(:,1,i_tau)+f_tau(:,2,i_tau)
	    bt(:,i_tau)=f_tau(:,1,i_tau)-f_tau(:,2,i_tau)
	  enddo
      call dgemm('n','n',2*n,n_nu+1,n_tau/2+1,1.d0,at,2*n,
     &           ab_nu_from_tau(0,0,1),n_tau/2+1+dp_tau_to_nu,0.d0,an,
     &           2*n)
      call dgemm('n','n',2*n,n_nu+1,n_tau/2+1,1.d0,bt,2*n,
     &           ab_nu_from_tau(0,0,2),n_tau/2+1+dp_tau_to_nu,0.d0,bn,
     &           2*n)
	  do i=1,n
        do i_nu=0,n_nu
	      f_nu(i,i_nu,1)=an(i,i_nu)-(0.d0,1.d0)*bn(i,i_nu)
	      f_nu(i,i_nu,2)=an(i,i_nu)+(0.d0,1.d0)*bn(i,i_nu)
	    enddo
	  enddo
	  deallocate(at,bt,an,bn)
      end




      subroutine tau_to_nu_ba3(f_tn,m0,n0,m,n,key)
c   ------- Transform TAU ---> NU representation ----------------------
c	  F_tau - is stored in A form
c     key = 1 for COS transform
c     key = 2 for SIN transform      
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,m,m0,key
	  real*8, intent(inout) :: f_tn(m0,n0,ndim3_tn)
      integer :: ind_tau,ind_nu,j,i_nu
	  real*8, allocatable :: at(:,:),an(:,:)
	  allocate(at(m,ndim3_tau))
	  allocate(an(m,0:n_nu))
      do j=1,n
        do ind_tau=1,ndim3_tau
	      at(:,ind_tau)=f_tn(1:m,j,ind_tau)
	    enddo
        call dgemm('n','n',m,n_nu+1,ndim3_tau,1.d0,at,m,
     &             ab_nu_from_tau(me_t*ndim3_tau,0,key),
     &             n_tau/2+1+dp_tau_to_nu,0.d0,an,m)
        if(nproc_t/=1) call dgop(an,m*(n_nu+1),'  +',comm_t)
        do ind_nu=1,ndim3_nu
	      i_nu=me_t*ndim3_nu+ind_nu-1
	      f_tn(1:m,j,ind_nu)=an(:,i_nu)
	    enddo
      enddo
	  deallocate(at,an)
      end




      subroutine tau_to_nu_ba3ut(f_tn,n0,n,key)
c   ------- Transform TAU ---> NU representation ----------------------
c	  F_tau - is stored in A form
c     key = 1 for COS transform
c     key = 2 for SIN transform
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,key
	  real*8, intent(inout) :: f_tn(n0,ndim3_tn)
      integer :: ind_nu,i_nu,nblock,nstep,i1,i2,istep
	  real*8, allocatable :: an(:,:)
	  nblock=1024
	  nstep=n/nblock
	  if(nstep*nblock<n) nstep=nstep+1
	  allocate(an(nblock,0:n_nu))
	  i1=1
	  i2=nblock
	  do istep=1,nstep
	    if(istep==nstep) i2=n-i1+1
        call dgemm('n','n',i2,n_nu+1,ndim3_tau,1.d0,f_tn(i1,1),n0,
     &             ab_nu_from_tau(me_t*ndim3_tau,0,key),
     &             n_tau/2+1+dp_tau_to_nu,0.d0,an,nblock)
        if(nproc_t/=1) call dgop(an,nblock*(n_nu+1),'  +',comm_t)
        do ind_nu=1,ndim3_nu
	      i_nu=me_t*ndim3_nu+ind_nu-1
	      f_tn(i1:i1+i2-1,ind_nu)=an(1:i2,i_nu)
        enddo
        i1=i1+nblock
      enddo
	  deallocate(an)
      end



      subroutine tau_to_nu_re(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation ----------------------
c   ---- We transform totally general real F_TAU with N elements ---
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [0 : n_nu] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_nu(n,ndim3_nu)
      integer :: ind_tau,ind_nu,i,i_nu
	real*8, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(bt(n,ndim3_tau))
	allocate(an(n,0:n_nu))
	allocate(bn(n,0:n_nu))
      do ind_tau=1,ndim3_tau
	  at(:,ind_tau)=f_tau(:,1,ind_tau)+f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	enddo
      call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,at,n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,1),
     &           n_tau/2+1+dp_tau_to_nu,0.d0,an,n)
      call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,bt,n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,2),
     &           n_tau/2+1+dp_tau_to_nu,0.d0,bn,n)
      if(nproc_t/=1) then
        call dgop(an,n*(n_nu+1),'  +',comm_t)
        call dgop(bn,n*(n_nu+1),'  +',comm_t)
      endif
	do i=1,n
	  do ind_nu=1,ndim3_nu
	    i_nu=me_t*ndim3_nu+ind_nu-1
	    f_nu(i,ind_nu)=dcmplx(an(i,i_nu),bn(i,i_nu))
	  enddo
	enddo
	deallocate(at,bt,an,bn)
      end



      subroutine tau_to_nu_aa_g(p_nu,p_tau,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  real*8, intent(in) :: p_tau(n,ndim3_tau)
	  real*8, intent(out) :: p_nu(n,ndim3_nu)
      integer :: ind_nu,i_nu
	  real*8, allocatable :: at(:,:),an(:,:)
	  allocate(at(n,ndim3_tau))
	  allocate(an(n,0:n_nu))
	  at=2.d0*p_tau
      an=0.d0
      call dgemm('n','n',n,n_nu+1,ndim3_tau,1.d0,at,n,
     &           ab_nu_from_tau(me_t*ndim3_tau,0,1),n_tau/2+1,0.d0,
     &           an,n)
	  if(nproc_t/=1) call dgop(an,n*(n_nu+1),'  +',comm_t)
	  do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
	    p_nu(:,ind_nu)=an(:,i_nu)
      enddo
      deallocate(at,an)
      end
