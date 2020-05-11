      subroutine from_nu_nu_to_tau_nu_g2_w2(n,ind_nu,nom,f_nu,f_tau)
c   ------- Transform T(v';v'-v) ----> AT(t;v) --------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,nom,ind_nu
      complex*16, intent(in) :: f_nu(n,nom)
      complex*16, intent(out) :: f_tau(n,0:n_tau)
      call zgemm('n','n',n,n_tau+1,nom,(1.d0,0.d0),f_nu,n,
     &           tau_nu_from_nu_nu_1(1,0,ind_nu),n_nu2_max,(0.d0,0.d0),
     &           f_tau,n)
      end
      
      
      
      
      subroutine from_nu_nu_to_tau_nu_re(n,ind_nu,nom,f_nu,f_tau)
c   ------- Transform T(v';v'-v) ----> AT(t;v) --------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,nom,ind_nu
      real*8, intent(in) :: f_nu(n,nom)
      complex*16, intent(out) :: f_tau(n,0:n_tau/2)
      real*8, allocatable :: a(:,:),b(:,:)
      allocate(a(nom,0:n_tau/2))
      allocate(b(n,0:n_tau/2))
      a=real(tau_nu_from_nu_nu_1(1:nom,0:n_tau/2,ind_nu))
      call dgemm('n','n',n,n_tau/2+1,nom,1.d0,f_nu,n,a,nom,0.d0,b,n)
      f_tau=b
      a=imag(tau_nu_from_nu_nu_1(1:nom,0:n_tau/2,ind_nu))
      call dgemm('n','n',n,n_tau/2+1,nom,1.d0,f_nu,n,a,nom,0.d0,b,n)
      f_tau=f_tau+(0.d0,1.d0)*b
      deallocate(a,b)
      end