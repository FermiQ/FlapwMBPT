      subroutine tau_to_omega_spl(f_om,f_tau,df_tau,n,key)
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
	  real*8, intent(in) :: f_tau(n,0:n_tau),df_tau(n,3)
	  complex*16, intent(out) :: f_om(n,0:n_omega)
      integer :: i_tau
      real*8, allocatable :: y(:,:),tmp(:,:)
      allocate(y(n,0:n_tau/2+key))
      allocate(tmp(n,0:n_omega))
c ----------------- A-part -------------------------------------
      do i_tau=0,n_tau/2
        y(:,i_tau)=f_tau(:,i_tau)-f_tau(:,n_tau-i_tau)
      enddo
      if(key==2) then
        y(:,n_tau/2+1)=df_tau(:,1)+df_tau(:,3)
        y(:,n_tau/2+2)=2.d0*df_tau(:,2)
      endif
      call dgemm('n','n',n,n_omega+1,n_tau/2+1+key,1.d0,y,n,
     &           ab_omega_from_tau(0,0,1),n_tau/2+1+key,0.d0,tmp,n)
      f_om=tmp
c ----------------- B-part -------------------------------------
      do i_tau=0,n_tau/2
        y(:,i_tau)=f_tau(:,i_tau)+f_tau(:,n_tau-i_tau)
      enddo
      if(key==2) then
        y(:,n_tau/2+1)=df_tau(:,1)-df_tau(:,3)
        y(:,n_tau/2+2)=0.d0
      endif
      call dgemm('n','n',n,n_omega+1,n_tau/2+1+key,1.d0,y,n,
     &           ab_omega_from_tau(0,0,2),n_tau/2+1+key,0.d0,tmp,n)
      f_om=f_om+(0.d0,1.d0)*tmp
      deallocate(y,tmp)
      end
      
      subroutine tau_to_omega_baa(f_omega,n0,f_tau,n1,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	  F_tau - hermitian, but is stored in square form = A;B
c     Output:
c     F_omega is complex*16 and stored in A,B form = (Re,Im)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f_tau(n1,n1,2,ndim3_tau)
	  real*8, intent(out) :: f_omega(n0,n0,2,ndim3_omega)
      integer :: nn,i,j,ij,ind_tau,ind_omega,it,i_omega
	  real*8, allocatable :: at(:,:),ao(:,:)
	  nn=n*n
	  allocate(at(nn,ndim3_tau))
	  allocate(ao(nn,0:n_omega))
c -------------------------- COS + SIN transformations ------------------
      do it=1,2
	    ao=0.d0
	    do ind_tau=1,ndim3_tau
	      ij=0
	      do j=1,n
            do i=1,n
	          ij=ij+1
	          at(ij,ind_tau)=f_tau(i,j,it,ind_tau)
	        enddo
          enddo
        enddo
        call dgemm('n','n',nn,n_omega+1,ndim3_tau,1.d0,at,nn,
     &             ab_omega_from_tau(me_t*ndim3_tau,0,it),
     &             n_tau/2+1+dsig_tau_to_omega,0.d0,ao,nn)
	    if(nproc_t/=1) call dgop(ao,nn*(n_omega+1),'  +',comm_t)
	    ij=0
	    do j=1,n
          do i=1,n
	        ij=ij+1
	        do ind_omega=1,ndim3_omega
	          i_omega=me_t*ndim3_omega+ind_omega-1
	          f_omega(i,j,it,ind_omega)=ao(ij,i_omega)
	        enddo
          enddo
        enddo
      enddo   !! over it
	  deallocate(at,ao)
      end



      subroutine tau_to_omega_vrt(f_om,f_tau,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_tau(n,2,ndim3_tau)
	  complex*16, intent(out) :: f_om(n,0:n_omega)
      integer :: ind_tau
	  complex*16, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
	  allocate(at(n,ndim3_tau))
	  allocate(ao(n,0:n_omega))
	  allocate(bt(n,ndim3_tau))
	  allocate(bo(n,0:n_omega))
	  do ind_tau=1,ndim3_tau
	    at(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	    bt(:,ind_tau)=(0.d0,1.d0)*(f_tau(:,1,ind_tau)
     &	                          +f_tau(:,2,ind_tau))
	  enddo
      call dgemm('n','n',2*n,n_omega+1,ndim3_tau,1.d0,at,2*n,
     &           ab_omega_from_tau(me_t*ndim3_tau,0,1),
     &           n_tau/2+dsig_tau_to_omega+1,0.d0,ao,2*n)
      call dgemm('n','n',2*n,n_omega+1,ndim3_tau,1.d0,bt,2*n,
     &           ab_omega_from_tau(me_t*ndim3_tau,0,2),
     &           n_tau/2+dsig_tau_to_omega+1,0.d0,bo,2*n)
      f_om=ao+bo
      deallocate(at,bt,ao,bo)
	  if(nproc_t/=1) call dgop(f_om,2*n*(n_omega+1),'  +',comm_t)
      end
      
      subroutine tau_to_omega_baac(f_omega,n0,f_tau,n1,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	  F_tau - hermitian, but is stored in square form = A;B
c     Output:
c     F_omega is complex*16 and stored in A,B form = (Re,Im)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f_tau(n1,n1,2,ndim3_tau)
	  real*8, intent(out) :: f_omega(n0,n0,2,ndim3_omega)
      integer :: nn,i,j,ij,ind_tau,ind_omega,it,i_omega
	  real*8, allocatable :: at(:,:),ao(:,:)
	  nn=n*n
	  allocate(at(nn,ndim3_tau))
	  allocate(ao(nn,0:n_omega))
c -------------------------- COS + SIN transformations ------------------
      do it=1,2
	    ao=0.d0
	    do ind_tau=1,ndim3_tau
	      ij=0
	      do j=1,n
            do i=1,n
	          ij=ij+1
	          at(ij,ind_tau)=f_tau(i,j,it,ind_tau)
	        enddo
          enddo
        enddo
        call dgemm('n','n',nn,n_omega+1,ndim3_tau,1.d0,at,nn,
     &             ab_omega_from_tau(me_t*ndim3_tau,0,it),
     &             n_tau/2+1+dsig_tau_to_omega,0.d0,ao,nn)
	  if(nproc_t/=1) call dgop(ao,nn*(n_omega+1),'  +',comm_t)
	  ij=0
	  do j=1,n
          do i=1,n
	        ij=ij+1
	        do ind_omega=1,ndim3_omega
              i_omega=me_t*ndim3_omega+ind_omega-1
	          f_omega(i,j,it,ind_omega)=ao(ij,i_omega)
	        enddo
          enddo
        enddo
      enddo   !! over it
	  deallocate(at,ao)
      end
      
      subroutine tau_to_omega_seq(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	  F_tau - C*16, N elements in (0:n_tau) form
c     Output:
c     F_omega - C*16, N elemnents (-+ omega)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_tau(n,0:n_tau)
	  complex*16, intent(out) :: f_omega(n,0:n_omega,2)
      integer :: i_tau,i_omega,i1_tau
	  complex*16, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
	  allocate(at(n,0:n_tau/2))
	  allocate(ao(n,0:n_omega))
	  allocate(bt(n,0:n_tau/2))
	  allocate(bo(n,0:n_omega))
	  do i_tau=0,n_tau/2
	    i1_tau=n_tau-i_tau
	    at(:,i_tau)=f_tau(:,i_tau)-f_tau(:,i1_tau)
	    bt(:,i_tau)=(0.d0,1.d0)*(f_tau(:,i_tau)+f_tau(:,i1_tau))
	  enddo
	  call dgemm('n','n',2*n,n_omega+1,n_tau/2+1,1.d0,at,2*n,
     &           ab_omega_from_tau(0,0,1),n_tau/2+1,0.d0,ao,2*n)
	  call dgemm('n','n',2*n,n_omega+1,n_tau/2+1,1.d0,bt,2*n,
     &           ab_omega_from_tau(0,0,2),n_tau/2+1,0.d0,bo,2*n)
	  do i_omega=0,n_omega
	    f_omega(:,i_omega,1)=ao(:,i_omega)-bo(:,i_omega)
	    f_omega(:,i_omega,2)=ao(:,i_omega)+bo(:,i_omega)
	  enddo
	  deallocate(at,ao,bt,bo)
      end




      subroutine tau_to_omega_abb(f_om,n0,f_tau,n1,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,n0,n1
	  real*8, intent(in) :: f_tau(n1,n1,2,ndim3_tau)
	  complex*16, intent(out) :: f_om(n0,n0,ndim3_omega)
      integer :: ind_tau,nn,ind,ind_omega,i_omega
      real*8, allocatable :: y(:,:,:),tmp(:,:,:)
      allocate(y(n,n,ndim3_tau))
      allocate(tmp(n,n,0:n_omega))
c ----------------- A-part -------------------------------------
      do ind_tau=1,ndim3_tau
        y(:,:,ind_tau)=f_tau(1:n,1:n,1,ind_tau)
     &                -f_tau(1:n,1:n,2,ind_tau)
      enddo
      nn=n*n
      ind=me_t*ndim3_tau
      call dgemm('n','n',nn,n_omega+1,ndim3_tau,1.d0,y,nn,
     &           ab_omega_from_tau(ind,0,1),n_tau/2+1+dsig_tau_to_omega,
     &           0.d0,tmp,nn)
      if(nproc_t/=1) call dgop(tmp,nn*(n_omega+1),'  +',comm_t)
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        f_om(1:n,1:n,ind_omega)=tmp(:,:,i_omega)
      enddo
c ----------------- B-part -------------------------------------
      do ind_tau=1,ndim3_tau
        y(:,:,ind_tau)=f_tau(1:n,1:n,1,ind_tau)
     &                +f_tau(1:n,1:n,2,ind_tau)
      enddo
      call dgemm('n','n',nn,n_omega+1,ndim3_tau,1.d0,y,nn,
     &           ab_omega_from_tau(ind,0,2),n_tau/2+1+dsig_tau_to_omega,
     &           0.d0,tmp,nn)
      if(nproc_t/=1) call dgop(tmp,nn*(n_omega+1),'  +',comm_t)
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        f_om(1:n,1:n,ind_omega)=f_om(1:n,1:n,ind_omega)
     &                         +(0.d0,1.d0)*tmp(:,:,i_omega)
      enddo
      deallocate(y,tmp)
      end




      subroutine tau_to_omega_abbc(f_om,n0,f_tau,n1,n)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n0,n1
      complex*16, intent(in) :: f_tau(n1,n1,2,ndim3_tau)
      complex*16, intent(out) :: f_om(n0,n0,ndim3_omega)
      integer :: ind_tau,nn,ind,ind_omega,i_omega
      complex*16, allocatable :: y(:,:,:),tmp(:,:,:)
      allocate(y(n,n,ndim3_tau))
      allocate(tmp(n,n,0:n_omega))
c ----------------- A-part -------------------------------------
      do ind_tau=1,ndim3_tau
        y(:,:,ind_tau)=f_tau(1:n,1:n,1,ind_tau)
     &                -f_tau(1:n,1:n,2,ind_tau)
      enddo
      nn=2*n*n
      ind=me_t*ndim3_tau
      call dgemm('n','n',nn,n_omega+1,ndim3_tau,1.d0,y,nn,
     &           ab_omega_from_tau(ind,0,1),n_tau/2+1+dsig_tau_to_omega,
     &           0.d0,tmp,nn)
      if(nproc_t/=1) call dgop(tmp,nn*(n_omega+1),'  +',comm_t)
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        f_om(1:n,1:n,ind_omega)=tmp(:,:,i_omega)
      enddo
c ----------------- B-part -------------------------------------
      do ind_tau=1,ndim3_tau
        y(:,:,ind_tau)=f_tau(1:n,1:n,1,ind_tau)
     &                +f_tau(1:n,1:n,2,ind_tau)
      enddo
      call dgemm('n','n',nn,n_omega+1,ndim3_tau,1.d0,y,nn,
     &           ab_omega_from_tau(ind,0,2),n_tau/2+1+dsig_tau_to_omega,
     &           0.d0,tmp,nn)
      if(nproc_t/=1) call dgop(tmp,nn*(n_omega+1),'  +',comm_t)
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        f_om(1:n,1:n,ind_omega)=f_om(1:n,1:n,ind_omega)
     &                         +(0.d0,1.d0)*tmp(:,:,i_omega)
      enddo
      deallocate(y,tmp)
      end
      
      subroutine tau_to_omega_kk(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	  F_tau - C*16, N elements in (0:n_tau) form
c     Output:
c     F_omega - C*16, N elemnents (-+ omega)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n
	  complex*16, intent(in) :: f_tau(n,2,0:n_tau/2)
	  complex*16, intent(out) :: f_omega(n,0:n_omega,2)
      integer :: i_omega,i_tau
	  complex*16, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
	  allocate(at(n,0:n_tau/2))
	  allocate(ao(n,0:n_omega))
	  allocate(bt(n,0:n_tau/2))
	  allocate(bo(n,0:n_omega))
	  do i_tau=0,n_tau/2
	    at(:,i_tau)=f_tau(:,1,i_tau)-f_tau(:,2,i_tau)
	    bt(:,i_tau)=(0.d0,1.d0)*(f_tau(:,1,i_tau)+f_tau(:,2,i_tau))
	  enddo
	  call dgemm('n','n',2*n,n_omega+1,n_tau/2+1,1.d0,at,2*n,
     &           ab_omega_from_tau(0,0,1),n_tau/2+1,0.d0,ao,2*n)
	  call dgemm('n','n',2*n,n_omega+1,n_tau/2+1,1.d0,bt,2*n,
     &	       ab_omega_from_tau(0,0,2),n_tau/2+1,0.d0,bo,2*n)
	  do i_omega=0,n_omega
	    f_omega(:,i_omega,1)=ao(:,i_omega)-bo(:,i_omega)
	    f_omega(:,i_omega,2)=ao(:,i_omega)+bo(:,i_omega)
	  enddo
	  deallocate(at,ao,bt,bo)
      end


      subroutine tau_to_omega_c(f_om,f_tau,n)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      complex*16, intent(in) :: f_tau(n,0:n_tau)
      complex*16, intent(out) :: f_om(n,0:n_omega)
      integer :: i_tau
      complex*16, allocatable :: y(:,:),tmp(:,:)
      allocate(y(n,0:n_tau/2))
      allocate(tmp(n,0:n_omega))
c ----------------- A-part -------------------------------------
      do i_tau=0,n_tau/2
        y(:,i_tau)=f_tau(:,i_tau)-f_tau(:,n_tau-i_tau)
      enddo
      call dgemm('n','n',2*n,n_omega+1,n_tau/2+1,1.d0,y,2*n,
     &           ab_omega_from_tau(0,0,1),n_tau/2+1,0.d0,tmp,2*n)
      f_om=tmp
c ----------------- B-part -------------------------------------
      do i_tau=0,n_tau/2
        y(:,i_tau)=f_tau(:,i_tau)+f_tau(:,n_tau-i_tau)
      enddo
      call dgemm('n','n',2*n,n_omega+1,n_tau/2+1,1.d0,y,2*n,
     &           ab_omega_from_tau(0,0,2),n_tau/2+1,0.d0,tmp,2*n)
      f_om=f_om+(0.d0,1.d0)*tmp
      deallocate(y,tmp)
      end

      subroutine tau_to_omega_ba1(f_omega,n0,f_tau,n1,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c     F_tau - hermitian, but is stored in square form = A;B
c     Output:
c     F_omega is complex*16 and stored in A,B form = (Re,Im)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n0,n1
      real*8, intent(in) :: f_tau(n1,2,ndim3_tau)
      real*8, intent(out) :: f_omega(n0,2,ndim3_omega)
      integer :: i,ij,ind_tau,ind_omega,it,i_omega
      real*8, allocatable :: at(:,:),ao(:,:)
      allocate(at(n,ndim3_tau))
      allocate(ao(n,0:n_omega))
c -------------------------- COS + SIN transformations ------------------
      do it=1,2
        ao=0.d0
        do ind_tau=1,ndim3_tau
          do i=1,n
            at(i,ind_tau)=f_tau(i,it,ind_tau)
          enddo
        enddo
        call dgemm('n','n',n,n_omega+1,ndim3_tau,1.d0,at,n,
     &             ab_omega_from_tau(me_t*ndim3_tau,0,it),
     &             n_tau/2+1+dsig_tau_to_omega,0.d0,ao,n)
        if(nproc_t/=1) call dgop(ao,n*(n_omega+1),'  +',comm_t)
        ij=0
        do i=1,n
          do ind_omega=1,ndim3_omega
            i_omega=me_t*ndim3_omega+ind_omega-1
            f_omega(i,it,ind_omega)=ao(i,i_omega)
          enddo
        enddo
      enddo   !! over it
      deallocate(at,ao)
      end
