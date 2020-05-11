      subroutine nu_to_tau_spl(p_nu,p_tau,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  real*8, intent(in) :: p_nu(n,0:n_nu)
	  real*8, intent(out) :: p_tau(n,0:n_tau/2)
      call dgemm('n','n',n,n_tau/2+1,n_nu+1,1.d0,p_nu,n,
     &           tau_from_nu(0,0,1),n_nu+1,0.d0,p_tau,n)
      end
      
      
      subroutine nu_to_tau_c_spl(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is C*16 nonsymmetric; stored in N elements [0:n_nu]
c     Output:
c	  F_tau - R*8 nonsymmetric; stored in N elements (t;beta-t)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_nu(n,ndim3_nu)
	  real*8, intent(out) :: f_tau(n,2,ndim3_tau)
      integer :: ind_nu,ind_tau,i_nu
	  real*8, allocatable :: at(:,:),an(:,:),bt(:,:),bn(:,:)
	  allocate(at(n,ndim3_tau))
	  allocate(an(n,0:n_nu))
	  allocate(bt(n,ndim3_tau))
	  allocate(bn(n,0:n_nu))
	  an=0.d0
	  bn=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
	    an(:,i_nu)=real(f_nu(:,ind_nu))
	    bn(:,i_nu)=imag(f_nu(:,ind_nu))
      enddo
	  if(nproc_t/=1) then
	    call dgop(an,n*(n_nu+1),'  +',comm_t)
	    call dgop(bn,n*(n_nu+1),'  +',comm_t)
	  endif
	  call dgemm('n','n',n,ndim3_tau,n_nu+1,1.d0,
     &           an(1,0),n,tau_from_nu(0,me_t*ndim3_tau,1),n_nu+1,
     &           0.d0,at,n)
	  call dgemm('n','n',n,ndim3_tau,n_nu+1,1.d0,
     &           bn(1,0),n,tau_from_nu(0,me_t*ndim3_tau,2),n_nu+1,
     &           0.d0,bt,n)
	  do ind_tau=1,ndim3_tau
	    f_tau(:,1,ind_tau)=0.5d0*(at(:,ind_tau)+bt(:,ind_tau))
	    f_tau(:,2,ind_tau)=0.5d0*(at(:,ind_tau)-bt(:,ind_tau))
	  enddo
	  deallocate(an,at,bn,bt)
      end




      subroutine nu_to_tau_ba(p_nu,n0,p_tau,n1,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: p_nu(n0,n0,ndim3_nu)
	  real*8, intent(out) :: p_tau(n1,n1,ndim3_tau)
      integer :: j,ind_nu,ind_tau,i_tau
	  real*8, allocatable :: at(:,:),an(:,:)
	  allocate(at(n,0:n_tau/2))
	  allocate(an(n,ndim3_nu))
      do j=1,n
        do ind_nu=1,ndim3_nu
	      an(:,ind_nu)=p_nu(1:n,j,ind_nu)
	    enddo
        call dgemm('n','n',n,n_tau/2+1,ndim3_nu,1.d0,an,n,
     &             tau_from_nu(me_t*ndim3_nu,0,1),n_nu+1,0.d0,at,n)
        if(nproc_t/=1) call dgop(at,n*(n_tau/2+1),'  +',comm_t)
        do ind_tau=1,ndim3_tau
	      i_tau=me_t*ndim3_tau+ind_tau-1
	      p_tau(1:n,j,ind_tau)=at(:,i_tau)
	    enddo
      enddo
      deallocate(at,an)
      end




      subroutine nua_to_tau_ba(p_nu,n0,p_tau,n1,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: p_nu(n0,n0,0:n_nu)
	  real*8, intent(out) :: p_tau(n1,n1,ndim3_tau)
      integer :: j,ind_tau,i_nu
	  real*8, allocatable :: at(:,:),an(:,:)
	  allocate(at(n,ndim3_tau))
	  allocate(an(n,0:n_nu))
      do j=1,n
        do i_nu=0,n_nu
	      an(:,i_nu)=p_nu(1:n,j,i_nu)
	    enddo
        call dgemm('n','n',n,ndim3_tau,n_nu+1,1.d0,an,n,
     &             tau_from_nu(0,me_t*ndim3_tau,1),n_nu+1,
     &             0.d0,at,n)
        do ind_tau=1,ndim3_tau
	      p_tau(1:n,j,ind_tau)=at(:,ind_tau)
	    enddo
      enddo
      end
      
      
      subroutine nu_to_tau_cc(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is given by its F(-v) and F(v) components (N elements)
c     Output:
c	  F_tau - C*16 nonsymmetric; stored in N elements (t;beta-t)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_nu(n,ndim3_nu,2)
	  complex*16, intent(out) :: f_tau(n,2,ndim3_tau)
      integer :: ind_nu,ind_tau,i_nu
	  complex*16, allocatable :: at(:,:),an(:,:),bt(:,:),bn(:,:)
	  allocate(at(n,ndim3_tau))
	  allocate(an(n,0:n_nu))
	  allocate(bt(n,ndim3_tau))
	allocate(bn(n,0:n_nu))
	an=(0.d0,0.d0)
	bn=(0.d0,0.d0)
	do ind_nu=1,ndim3_nu
	  i_nu=me_t*ndim3_nu+ind_nu-1
	  an(:,i_nu)=0.5d0*(f_nu(:,ind_nu,1)+f_nu(:,ind_nu,2))
	  bn(:,i_nu)=(0.d0,0.5d0)*(f_nu(:,ind_nu,1)-f_nu(:,ind_nu,2))
	enddo
	if(nproc_t/=1) then
	  call dgop(an,2*n*(n_nu+1),'  +',comm_t)
	  call dgop(bn,2*n*(n_nu+1),'  +',comm_t)
	endif
	call dgemm('n','n',2*n,ndim3_tau,n_nu+1,1.d0,
     &           an(1,0),2*n,tau_from_nu(0,me_t*ndim3_tau,1),n_nu+1,
     &           0.d0,at,2*n)
	call dgemm('n','n',2*n,ndim3_tau,n_nu+1,1.d0,
     &           bn(1,0),2*n,tau_from_nu(0,me_t*ndim3_tau,2),n_nu+1,
     &           0.d0,bt,2*n)
	do ind_tau=1,ndim3_tau
	  f_tau(:,1,ind_tau)=0.5d0*(at(:,ind_tau)+bt(:,ind_tau))
	  f_tau(:,2,ind_tau)=0.5d0*(at(:,ind_tau)-bt(:,ind_tau))
	enddo
	deallocate(at,an,bt,bn)
      end
      
      
      subroutine nu_to_tau_seq(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is given by its F(-v) and F(v) components (N elements)
c     Output:
c	F_tau - C*16 nonsymmetric; stored in N elements (t;beta-t)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_nu(n,0:n_nu,2)
	complex*16, intent(out) :: f_tau(n,0:n_tau)
      integer :: i_tau,i_nu,i1_tau
	complex*16, allocatable :: at(:,:),an(:,:),bt(:,:),bn(:,:)
	allocate(at(n,0:n_tau/2))
	allocate(an(n,0:n_nu))
	allocate(bt(n,0:n_tau/2))
	allocate(bn(n,0:n_nu))
	an=(0.d0,0.d0)
	bn=(0.d0,0.d0)
	do i_nu=0,n_nu
	  an(:,i_nu)=0.5d0*(f_nu(:,i_nu,1)+f_nu(:,i_nu,2))
	  bn(:,i_nu)=(0.d0,0.5d0)*(f_nu(:,i_nu,1)-f_nu(:,i_nu,2))
	enddo
	call dgemm('n','n',2*n,n_tau/2+1,n_nu+1,1.d0,an,2*n,
     &           tau_from_nu(0,0,1),n_nu+1,0.d0,at,2*n)
	call dgemm('n','n',2*n,n_tau/2+1,n_nu+1,1.d0,bn,2*n,
     &           tau_from_nu(0,0,2),n_nu+1,0.d0,bt,2*n)
	do i_tau=0,n_tau/2
	  i1_tau=n_tau-i_tau
	  f_tau(:,i_tau)=0.5d0*(at(:,i_tau)+bt(:,i_tau))
	  f_tau(:,i1_tau)=0.5d0*(at(:,i_tau)-bt(:,i_tau))
	enddo
	deallocate(at,an,bt,bn)
      end
      
      
      subroutine nu_to_tau_half(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c   --- Sums up only v>=0 components ----------------------------
c     Input:
c     F_nu is given by its F(v) components (N elements)
c     Output:
c	F_tau - C*16 nonsymmetric; stored in N elements (t;beta-t)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_nu(n,ndim3_nu)
	complex*16, intent(out) :: f_tau(n,2,ndim3_tau)
      integer :: ind_nu,ind_tau,i_nu
	complex*16, allocatable :: at(:,:),an(:,:),bt(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(an(n,0:n_nu))
	allocate(bt(n,ndim3_tau))
	allocate(bn(n,0:n_nu))
	an=(0.d0,0.d0)
	bn=(0.d0,0.d0)
	do ind_nu=1,ndim3_nu
	  i_nu=me_t*ndim3_nu+ind_nu-1
	  an(:,i_nu)=0.5d0*f_nu(:,ind_nu)
	  bn(:,i_nu)=(0.d0,-0.5d0)*f_nu(:,ind_nu)
	enddo
	if(nproc_t/=1) then
	  call dgop(an,2*n*(n_nu+1),'  +',comm_t)
	  call dgop(bn,2*n*(n_nu+1),'  +',comm_t)
	endif
	call dgemm('n','n',2*n,ndim3_tau,n_nu+1,1.d0,
     &           an(1,0),2*n,tau_from_nu(0,me_t*ndim3_tau,1),n_nu+1,
     &           0.d0,at,2*n)
	call dgemm('n','n',2*n,ndim3_tau,n_nu+1,1.d0,
     &           bn(1,0),2*n,tau_from_nu(0,me_t*ndim3_tau,2),n_nu+1,
     &           0.d0,bt,2*n)
	do ind_tau=1,ndim3_tau
	  f_tau(:,1,ind_tau)=0.5d0*(at(:,ind_tau)+bt(:,ind_tau))
	  f_tau(:,2,ind_tau)=0.5d0*(at(:,ind_tau)-bt(:,ind_tau))
	enddo
	deallocate(at,an,bt,bn)
      end
      
      
      subroutine nu_to_tau_kk4(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is given by its F(-v) and F(v) components (N elements)
c     Output:
c	F_tau - C*16 nonsymmetric; stored in N elements (t;beta-t)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_nu(n,0:n_nu,2)
	complex*16, intent(out) :: f_tau(n,2,0:n_tau/2)
      integer :: i_tau,i_nu
	complex*16, allocatable :: at(:,:),an(:,:),bt(:,:),bn(:,:)
	allocate(at(n,0:n_tau/2))
	allocate(an(n,0:n_nu))
	allocate(bt(n,0:n_tau/2))
	allocate(bn(n,0:n_nu))
	an=(0.d0,0.d0)
	bn=(0.d0,0.d0)
	do i_nu=0,n_nu
	  an(:,i_nu)=0.5d0*(f_nu(:,i_nu,1)+f_nu(:,i_nu,2))
	  bn(:,i_nu)=(0.d0,0.5d0)*(f_nu(:,i_nu,1)-f_nu(:,i_nu,2))
	enddo
	call dgemm('n','n',2*n,n_tau/2+1,n_nu+1,1.d0,an(1,0),2*n,
     &           tau_from_nu(0,0,1),n_nu+1,0.d0,at,2*n)
	call dgemm('n','n',2*n,n_tau/2+1,n_nu+1,1.d0,bn(1,0),2*n,
     &           tau_from_nu(0,0,2),n_nu+1,0.d0,bt,2*n)
	do i_tau=0,n_tau/2
	  f_tau(:,1,i_tau)=0.5d0*(at(:,i_tau)+bt(:,i_tau))
	  f_tau(:,2,i_tau)=0.5d0*(at(:,i_tau)-bt(:,i_tau))
	enddo
	deallocate(at,an,bt,bn)
      end




      subroutine nu_to_tau_ba3(p_tn,m0,n0,m,n)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n0,m,m0
      real*8, intent(inout) :: p_tn(m0,n0,ndim3_tn)
      integer :: j,ind_nu,ind_tau,i_tau
      real*8, allocatable :: at(:,:),an(:,:)
      allocate(at(m,0:n_tau/2))
      allocate(an(m,ndim3_nu))
      do j=1,n
        do ind_nu=1,ndim3_nu
	      an(:,ind_nu)=p_tn(1:m,j,ind_nu)
        enddo
        call dgemm('n','n',m,n_tau/2+1,ndim3_nu,1.d0,an,m,
     &             tau_from_nu(me_t*ndim3_nu,0,1),n_nu+1,0.d0,at,m)
        if(nproc_t/=1) call dgop(at,m*(n_tau/2+1),'  +',comm_t)
        do ind_tau=1,ndim3_tau
          i_tau=me_t*ndim3_tau+ind_tau-1
	      p_tn(1:m,j,ind_tau)=at(:,i_tau)
	    enddo
      enddo
      deallocate(at,an)
      end




      subroutine nu_to_tau_ba3ut(p_tn,n0,n)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n0
      real*8, intent(inout) :: p_tn(n0,ndim3_tn)
      integer :: ind_tau,i_tau,nblock,nstep,istep,i1,i2
      real*8, allocatable :: at(:,:)
      nblock=1024
      nstep=n/nblock
      if(nstep*nblock<n) nstep=nstep+1
      allocate(at(nblock,0:n_tau/2))
      i1=1
      i2=nblock
      do istep=1,nstep
        if(istep==nstep) i2=n-i1+1
        call dgemm('n','n',i2,n_tau/2+1,ndim3_nu,1.d0,p_tn(i1,1),n0,
     &             tau_from_nu(me_t*ndim3_nu,0,1),n_nu+1,0.d0,at,nblock)
        if(nproc_t/=1) call dgop(at,nblock*(n_tau/2+1),'  +',comm_t)
        do ind_tau=1,ndim3_tau
	      i_tau=me_t*ndim3_tau+ind_tau-1
	      p_tn(i1:i1+i2-1,ind_tau)=at(1:i2,i_tau)
        enddo
        i1=i1+nblock
      enddo
      deallocate(at)
      end
      
       subroutine nu_to_tau_aa_r(f_nu,n0,f_tau,n1,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is real symmetric but stored in square real*8 form
c     Output:
c	  F_tau - real symmetric, but is stored in square form
c     F(tau) = F(Beta-tau)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f_nu(n0,n0,0:n_nu)
	  real*8, intent(out) :: f_tau(n1,n1,0:n_tau/2)
      integer :: nn,i,j,ij,i_nu,i_tau
	  real*8, allocatable :: at(:,:),an(:,:)
	  nn=n*(n+1)/2
	  allocate(at(nn,0:n_tau/2))
	  allocate(an(nn,0:n_nu))
	  do i_nu=0,n_nu
	    ij=0
	    do j=1,n
	      do i=1,j
	        ij=ij+1
	        an(ij,i_nu)=f_nu(i,j,i_nu)
	      enddo
	    enddo
	  enddo
	  call dgemm('n','n',nn,n_tau/2+1,n_nu+1,1.d0,an(1,0),nn,
     &           tau_from_nu(0,0,1),n_nu+1,0.d0,at,nn)
	  f_tau=0.d0
	  do i_tau=0,n_tau/2
	    ij=0
	    do j=1,n
	      do i=1,j
	        ij=ij+1
	        f_tau(i,j,i_tau)=0.5d0*at(ij,i_tau)
	        f_tau(j,i,i_tau)=f_tau(i,j,i_tau)
	      enddo
	    enddo
	  enddo
	  deallocate(at,an)
      end



      subroutine nu_to_tau_aa_g(p_nu,p_tau,n)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: p_nu(n,ndim3_nu)
	real*8, intent(out) :: p_tau(n,ndim3_tau)
      integer :: ind_tau,i_tau
	real*8, allocatable :: at(:,:)
	allocate(at(n,0:n_tau/2))
      at=0.d0
      call dgemm('n','n',n,n_tau/2+1,ndim3_nu,1.d0,p_nu,n,
     &           tau_from_nu(me_t*ndim3_nu,0,1),n_nu+1,0.d0,at,n)
	if(nproc_t/=1) call dgop(at,n*(n_tau/2+1),'  +',comm_t)
	do ind_tau=1,ndim3_tau
        i_tau=me_t*ndim3_tau+ind_tau-1
	  p_tau(:,ind_tau)=0.5d0*at(:,i_tau)
      enddo
      deallocate(at)
      end


      subroutine from_nu_to_tau0(n,f_nu,f_tau)
c   ------- Transform F(NU) ---> F(TAU=0)  --------------------------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: f_nu(n,0:n_nu)
      real*8, intent(out) :: f_tau(n)
      call dgemm('n','n',n,1,n_nu+1,1.d0,f_nu(1,0),n,
     &           tau_from_nu(0,0,1),n_nu+1,0.d0,f_tau,n)
      f_tau=f_tau/2.d0
      end
