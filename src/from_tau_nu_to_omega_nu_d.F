      subroutine from_tau_nu_to_omega_nu_d(n,ind_nu,nom,f_om,f_tau,key)
c   ------- Transform F(t;v) ---> F(-w+v;-w) : KEY = 1 -------
c   ------- Transform F(t;v) ---> F(w;w-v)  : KEY = 2 ------
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,key,ind_nu,nom
	  complex*16, intent(in) :: f_tau(n,0:n_tau,2)
	  complex*16, intent(out) :: f_om(n,nom)
	  integer :: i_tau,i1_tau
	  complex*16, allocatable :: c(:,:),a(:,:,:),b(:,:,:),d(:,:)
	  allocate(c(n,nom))
	  allocate(d(n,nom))
	  allocate(a(n,0:n_tau/2,2))
	  allocate(b(n,0:n_tau/2,2))
	  do i_tau=0,n_tau/2
	    i1_tau=n_tau-i_tau
	    a(:,i_tau,1)=f_tau(:,i_tau,1)-f_tau(:,i1_tau,1)
	    a(:,i_tau,2)=f_tau(:,i_tau,2)-f_tau(:,i1_tau,2)
	    b(:,i_tau,1)=f_tau(:,i_tau,1)+f_tau(:,i1_tau,1)
	    b(:,i_tau,2)=f_tau(:,i_tau,2)+f_tau(:,i1_tau,2)
	  enddo
	  if(key==2) then   ! for w >= v/2
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,a(1,0,1),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,1,ind_nu),n_tau/2+1,
     &             0.d0,d,2*n)
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,b(1,0,1),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,2,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        d=d+(0.d0,1.d0)*c
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,a(1,0,2),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,3,ind_nu),n_tau/2+1,
     &             1.d0,d,2*n)
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,b(1,0,2),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,4,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        d=d+(0.d0,1.d0)*c
	  else if(key==1) then   ! for w <= v/2
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,a(1,0,1),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,3,ind_nu),n_tau/2+1,
     &             0.d0,d,2*n)
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,b(1,0,1),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,4,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        d=d+(0.d0,-1.d0)*c
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,a(1,0,2),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,1,ind_nu),n_tau/2+1,
     &             1.d0,d,2*n)
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,b(1,0,2),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,2,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        d=d+(0.d0,-1.d0)*c
      endif
      f_om=d
      deallocate(c,a,b,d)
      end




      subroutine from_tau_nu_to_omega_nu_kk(n,ind_nu,nom,f_om,f_tau,key)
c   ------- Transform F(t;v) ---> F(-w+v;-w) : KEY = 1 -------
c   ------- Transform F(t;v) ---> F(w;w-v)  : KEY = 2 ------
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: n,key,ind_nu,nom
	  complex*16, intent(in) :: f_tau(n,2,0:n_tau/2,2)
	  complex*16, intent(out) :: f_om(n,nom)
	  integer :: i_tau
	  complex*16, allocatable :: c(:,:),a(:,:,:),b(:,:,:),d(:,:)
	  allocate(c(n,nom))
	  allocate(d(n,nom))
	  allocate(a(n,0:n_tau/2,2))
	  allocate(b(n,0:n_tau/2,2))
	  do i_tau=0,n_tau/2
	    a(:,i_tau,1)=f_tau(:,1,i_tau,1)-f_tau(:,2,i_tau,1)
	    a(:,i_tau,2)=f_tau(:,1,i_tau,2)-f_tau(:,2,i_tau,2)
	    b(:,i_tau,1)=f_tau(:,1,i_tau,1)+f_tau(:,2,i_tau,1)
	    b(:,i_tau,2)=f_tau(:,1,i_tau,2)+f_tau(:,2,i_tau,2)
	  enddo
	  if(key==2) then   ! for w >= v/2
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,a(1,0,1),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,1,ind_nu),n_tau/2+1,
     &             0.d0,d,2*n)
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,b(1,0,1),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,2,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        d=d+(0.d0,1.d0)*c
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,a(1,0,2),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,3,ind_nu),n_tau/2+1,
     &             1.d0,d,2*n)
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,b(1,0,2),2*n,
     &             omega_nu_from_tau_nu_ab(0,1,4,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        d=d+(0.d0,1.d0)*c
	  else if(key==1) then   ! for w <= v/2
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,a(1,0,1),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,3,ind_nu),n_tau/2+1,
     &             0.d0,d,2*n)
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,b(1,0,1),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,4,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        d=d+(0.d0,-1.d0)*c
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,a(1,0,2),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,1,ind_nu),n_tau/2+1,
     &             1.d0,d,2*n)
	    call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,b(1,0,2),2*n,
     &	           omega_nu_from_tau_nu_ab(0,1,2,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        d=d+(0.d0,-1.d0)*c
      endif
      f_om=d
      deallocate(c,a,b,d)
      end
