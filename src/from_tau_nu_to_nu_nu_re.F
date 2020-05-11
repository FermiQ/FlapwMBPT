      subroutine from_tau_nu_to_nu_nu_spl(n,ind_nu,nom,f_nu,f_tau)
c   ------- Transform AQ(t;v);BQ(t;v) ---> Q(v';v'-v) ---- N elements ----
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,ind_nu,nom
      real*8, intent(in) :: f_tau(n,0:n_tau/2,4)
      real*8, intent(out) :: f_nu(n,nom)
      real*8, allocatable :: a(:,:)
      allocate(a(n,nom))
c ---------------------- COS(v')t contribution ---------------------------
      call dgemm('n','n',n,nom,n_tau/2+1,1.d0,f_tau(1,0,1),n,
     &           nu_nu_from_tau_nu(0,1,1,1,ind_nu),n_tau/2+1,0.d0,a,
     &           n)
c ---------------------- SIN(v')t contribution ---------------------------
      call dgemm('n','n',n,nom,n_tau/2+1,1.d0,f_tau(1,0,2),n,
     &           nu_nu_from_tau_nu(0,1,2,1,ind_nu),n_tau/2+1,1.d0,a,
     &           n)
c ---------------------- COS(v'-v)t contribution ---------------------------
      call dgemm('n','n',n,nom,n_tau/2+1,1.d0,f_tau(1,0,3),n,
     &           nu_nu_from_tau_nu(0,1,1,2,ind_nu),n_tau/2+1,1.d0,a,
     &           n)
c ---------------------- SIN(v'-v)t contribution ---------------------------
      call dgemm('n','n',n,nom,n_tau/2+1,1.d0,f_tau(1,0,4),n,
     &           nu_nu_from_tau_nu(0,1,2,2,ind_nu),n_tau/2+1,1.d0,a,
     &           n)
      f_nu=a
      deallocate(a)
      end


      subroutine from_tau_nu_to_nu_nu_d1(n,tran,nom,f_nu,f_tau)
c   ------- Transform AQ(t;v);BQ(t;v) ---> Q(v';v'-v) ---- N elements ----
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n,nom
	real*8, intent(in) :: f_tau(n,0:n_tau/2,2),
     &                      tran(0:n_tau/2,nom,2)
	real*8, intent(out) :: f_nu(n,nom)
	real*8, allocatable :: a(:,:)
	allocate(a(n,nom))
c ---------------------- COS+COS_1 contribution ---------------------------
	call dgemm('n','n',n,nom,n_tau/2+1,1.d0,f_tau(1,0,1),n,tran,
     &	       n_tau/2+1,0.d0,a,n)
c ---------------------- SIN-SIN1 contribution ---------------------------
	call dgemm('n','n',n,nom,n_tau/2+1,1.d0,f_tau(1,0,2),n,
     &           tran(1,1,2),n_tau/2+1,1.d0,a,n)
      f_nu=a
      deallocate(a)
      end




      subroutine from_tau_nu_to_nu_nu_par(n,ind_nu,nd4,nm4,f_nu,f_tau)
c   ------- Transform AQ(t;v);BQ(t;v) ---> Q(v';v'-v) ---- N elements ----
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,ind_nu,nd4,nm4
      real*8, intent(in) :: f_tau(n,0:n_tau/2,4)
      real*8, intent(out) :: f_nu(n,nd4)
      real*8, allocatable :: a(:,:)
      allocate(a(n,nd4))
c ---------------------- COS(v')t contribution ---------------------------
      call dgemm('n','n',n,nd4,n_tau/2+1,1.d0,f_tau(1,0,1),n,
     &           nu_nu_from_tau_nu(0,nm4+1,1,1,ind_nu),n_tau/2+1,0.d0,a,
     &           n)
c ---------------------- SIN(v')t contribution ---------------------------
      call dgemm('n','n',n,nd4,n_tau/2+1,1.d0,f_tau(1,0,2),n,
     &           nu_nu_from_tau_nu(0,nm4+1,2,1,ind_nu),n_tau/2+1,1.d0,a,
     &           n)
c ---------------------- COS(v'-v)t contribution ---------------------------
      call dgemm('n','n',n,nd4,n_tau/2+1,1.d0,f_tau(1,0,3),n,
     &           nu_nu_from_tau_nu(0,nm4+1,1,2,ind_nu),n_tau/2+1,1.d0,a,
     &           n)
c ---------------------- SIN(v'-v)t contribution ---------------------------
      call dgemm('n','n',n,nd4,n_tau/2+1,1.d0,f_tau(1,0,4),n,
     &           nu_nu_from_tau_nu(0,nm4+1,2,2,ind_nu),n_tau/2+1,1.d0,a,
     &           n)
      f_nu=a
      deallocate(a)
      end