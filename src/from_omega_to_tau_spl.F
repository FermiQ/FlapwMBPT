      subroutine omega_to_tau(n,f_omega,f_tau,df_tau)
c   ------- Transform OMEGA ---> TAU representation ----------------
c     F(-w)=F*(w)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_omega(n,0:n_omega)
	  real*8, intent(out) :: f_tau(n,0:n_tau),df_tau(n,3)
      integer :: i_tau,i
	  real*8, allocatable :: at(:,:),bt(:,:),ao(:,:),da(:),db(:)
	  allocate(at(n,0:n_tau/2))
	  allocate(ao(n,0:n_omega))
	  allocate(bt(n,0:n_tau/2))
      allocate(da(n))
      allocate(db(n))
      ao=real(f_omega)
      call dgemm('n','n',n,n_tau/2+1,n_omega+1,1.d0,ao,n,
     &           tau_from_omega(0,0,1),n_omega+1,0.d0,at,n)
      do i=1,n
        da(i)=dot_product(ao(i,:),dtau_from_omega(:,1))
      enddo
      ao=imag(f_omega)
      call dgemm('n','n',n,n_tau/2+1,n_omega+1,1.d0,ao,n,
     &           tau_from_omega(0,0,2),n_omega+1,0.d0,bt,n)
      do i=1,n
        db(i)=dot_product(ao(i,:),dtau_from_omega(:,2))
      enddo
	  f_tau=0.d0
	  do i_tau=0,n_tau/2
	    f_tau(:,i_tau)=0.5d0*(at(:,i_tau)+bt(:,i_tau))
	    f_tau(:,n_tau-i_tau)=0.5d0*(bt(:,i_tau)-at(:,i_tau))
      enddo
      df_tau(:,1)=db/2
      df_tau(:,2)=da/2
      df_tau(:,3)=-df_tau(:,1)
      deallocate(ao,at,bt,da,db)
      end
      
      subroutine omega_to_tau_baa(f_omega,n0,f_tau,n1,n)
c   ------- Transform OMEGA ---> TAU  --------------------------
c     Input:
c     F_omega is complex*16 and stored in A,B form = (Re,Im)
c     Output:
c	F_tau - hermitian, but is stored in square form = A;B
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f_omega(n0,n0,2,ndim3_omega)
	  real*8, intent(out) :: f_tau(n1,n1,2,ndim3_tau)
      integer :: nn,i,j,ij,ind_tau,i_omega,ind_omega,it
	  real*8, allocatable :: at(:,:),ao(:,:)
	  nn=n*n
	  allocate(at(nn,ndim3_tau))
	  allocate(ao(nn,0:n_omega))
c -------------------------- COS + SIN transformations ------------------
      do it=1,2
	    ao=0.d0
	    do ind_omega=1,ndim3_omega
	      i_omega=me_t*ndim3_omega+ind_omega-1
	      ij=0
	      do j=1,n
	        do i=1,n
	          ij=ij+1
	          ao(ij,i_omega)=f_omega(i,j,it,ind_omega)
	        enddo
	      enddo
	    enddo
	    if(nproc_t/=1) call dgop(ao,nn*(n_omega+1),'  +',comm_t)
	    call dgemm('n','n',nn,ndim3_tau,n_omega+1,1.d0,ao,nn,
     &	           tau_from_omega(0,me_t*ndim3_tau,it),n_omega+1,
     &             0.d0,at,nn)
	    do ind_tau=1,ndim3_tau
	      ij=0
	      do j=1,n
	        do i=1,n
	          ij=ij+1
	          f_tau(i,j,it,ind_tau)=at(ij,ind_tau)
	        enddo
	      enddo
	    enddo
	  enddo   !! over it
	  deallocate(at,ao)
      end



      subroutine omega_to_tau_all(n,f_omega,f_tau)
c   ------- Transform OMEGA ---> TAU representation ----------------
c     F(-w)=F*(w)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_omega(n,ndim3_omega)
	  real*8, intent(out) :: f_tau(n,2,ndim3_tau)
      integer :: i_tau,m0,ind_tau
	  real*8, allocatable :: at(:,:),bt(:,:),ao(:,:)
	  allocate(at(n,0:n_tau/2))
	  allocate(ao(n,ndim3_omega))
	  allocate(bt(n,0:n_tau/2))
      m0=me_t*ndim3_omega
      at=0.d0
      bt=0.d0
      ao=real(f_omega)
      call dgemm('n','n',n,n_tau/2+1,ndim3_omega,1.d0,ao,n,
     &           tau_from_omega(m0,0,1),n_omega+1,0.d0,at,n)
      ao=imag(f_omega)
      call dgemm('n','n',n,n_tau/2+1,ndim3_omega,1.d0,ao,n,
     &           tau_from_omega(m0,0,2),n_omega+1,0.d0,bt,n)
	  if(nproc_t/=1) then
        call dgop(at,n*(n_tau/2+1),'  +',comm_t)
        call dgop(bt,n*(n_tau/2+1),'  +',comm_t)
      endif
	  do ind_tau=1,ndim3_tau
        i_tau=me_t*ndim3_tau+ind_tau-1
	    f_tau(:,1,ind_tau)=0.5d0*(at(:,i_tau)+bt(:,i_tau))
	    f_tau(:,2,ind_tau)=0.5d0*(bt(:,i_tau)-at(:,i_tau))
      enddo
      deallocate(ao,at,bt)
      end
