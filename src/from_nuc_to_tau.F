      subroutine from_nu_to_tau_sin(n,b_nu,b_t)
c   ------- Transform B(NU) ---> B(TAU)  --------------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: b_nu(n,0:n_nu)
	real*8, intent(out) :: b_t(n,ndim3_tau)
	call dgemm('n','n',n,ndim3_tau,n_nu+1,1.d0,
     &           b_nu(1,0),n,tau_from_nu(0,me_t*ndim3_tau,2),
     &           n_nu+1,0.d0,b_t,n)
      end
      
      
      subroutine from_nu_to_tau_cos(n,a_nu,a_tau)
c   ------- Transform A(NU) ---> A(TAU)  --------------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: a_nu(n,0:n_nu)
	real*8, intent(out) :: a_tau(n,ndim3_tau)
	call dgemm('n','n',n,ndim3_tau,n_nu+1,1.d0,
     &           a_nu(1,0),n,tau_from_nu(0,me_t*ndim3_tau,1),
     &           n_nu+1,0.d0,a_tau,n)
      end
      
      
      subroutine from_nu_to_tau_c(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is C*16 nonsymmetric; stored in N elements [0:n_nu]
c     Output:
c	F_tau - R*8 nonsymmetric; stored in N elements (t;beta-t)
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
	call from_nu_to_tau_cos(n,an,at)
	call from_nu_to_tau_sin(n,bn,bt)
	do ind_tau=1,ndim3_tau
	  f_tau(:,1,ind_tau)=0.5d0*(at(:,ind_tau)+bt(:,ind_tau))
	  f_tau(:,2,ind_tau)=0.5d0*(at(:,ind_tau)-bt(:,ind_tau))
	enddo
	deallocate(at,an,bt,bn)
      end