      subroutine phi_rpa_c(n,n0,ind_k,i_nu,p,v)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: n,n0,i_nu,ind_k
      complex*16, intent(in) :: p(n0,n0),v(n0,n0)
	  integer :: i,k
	  real*8 :: s
	  real*8, allocatable :: d(:)
	  complex*16, allocatable :: pc(:,:),vp(:,:),a(:,:)
	  allocate(pc(n,n))
	  allocate(vp(n,n))
	  allocate(a(n,n))
	  allocate(d(n))
	  k=n3_mpi_k(me_k+1)+ind_k
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),v,n0,p,n0,(0.d0,0.d0),vp,n)
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),p,n0,v,n0,(0.d0,0.d0),pc,n)
	  s=0.d0
	  do i=1,n
	    s=s+vp(i,i)
	  enddo
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),pc,n,vp,n,(0.d0,0.d0),a,n)
	  a=a-vp-pc
	  do i=1,n
	    a(i,i)=(1.d0,0.d0)+a(i,i)
	  enddo
	  call eig_val_only(n,n,a,d)
	  d=log(d)
	  s=s/2.d0+sum(d)/4.d0
	  phi_functional_c=phi_functional_c
     &                +et_sum_nu(i_nu)*s*wgt(k)*2.d0
	  deallocate(pc,vp,a,d)
      end
      
      
      subroutine heg_phi_rpa_c(p,v,s)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      real*8, intent(in) :: p,v
      real*8, intent(out) :: s
	  real*8 :: vp,a,d
	  vp=v*p
	  s=vp
	  a=(1.d0-vp)**2
	  d=log(a)
	  s=s/2.d0+d/4.d0
      end



      subroutine phi_rpa_c_r(n,n0,ind_k,i_nu,p,v)
c -   Evaluates correlation part of Luttinger Ward functional in RPA -----
      use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: n,n0,i_nu,ind_k
      real*8, intent(in) :: p(n0,n0),v(n0,n0)
	  integer :: i,k
	  real*8 :: s
	  real*8, allocatable :: d(:),pc(:,:),vp(:,:),a(:,:)
	  allocate(pc(n,n))
	  allocate(vp(n,n))
	  allocate(a(n,n))
	  allocate(d(n))
	  k=n3_mpi_k(me_k+1)+ind_k
	  call dgemm('n','n',n,n,n,1.d0,v,n0,p,n0,0.d0,vp,n)
	  call dgemm('n','n',n,n,n,1.d0,p,n0,v,n0,0.d0,pc,n)
	  s=0.d0
	  do i=1,n
	    s=s+vp(i,i)
	  enddo
	  call dgemm('n','n',n,n,n,1.d0,pc,n,vp,n,0.d0,a,n)
	  a=a-vp-pc
	  do i=1,n
	    a(i,i)=1.d0+a(i,i)
	  enddo
	  call eig_val_solver_r(n,n,a,d)
	  d=log(d)
	  s=s/2.d0+sum(d)/4.d0
	  phi_functional_c=phi_functional_c
     &                +et_sum_nu(i_nu)*s*wgt(k)*2.d0
	  deallocate(pc,vp,a,d)
      end
