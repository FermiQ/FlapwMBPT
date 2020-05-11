      subroutine from_omega_nu_to_tau_nu_ab(n,ind_nu,a_om,a_tau,key)
c   ------- Transform F(w;w-v) ---> A(t;v), B(t;v) --------
c     key=1   A
c     key=2   B
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,ind_nu,key
      complex*16, intent(in) :: a_om(n,n_omega1_max)
      complex*16, intent(out) :: a_tau(n,0:n_tau)
      integer :: nom,i_nu
	  i_nu=me_t*ndim3_nu+ind_nu-1
      nom=num_omega1_adapt(i_nu)
      call zgemm('n','n',n,n_tau+1,nom,(1.d0,0.d0),a_om,n,
     &	         tau_nu_from_omega_nu(1,0,ind_nu),n_omega1_max,
     &           (0.d0,0.d0),a_tau,n)
      if(key==1) a_tau=conjg(a_tau)
      end

      subroutine from_omega_nu_to_tau_nu(n,i_nu,a_om,a_tau,key)
c   ------- Transform F(w;w-v) ---> A(t;v), B(t;v) --------
c     key=1   A
c     key=2   B
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,i_nu,key
      complex*16, intent(in) :: a_om(n,n_omega1_max)
      complex*16, intent(out) :: a_tau(n,2,ndim3_tau)
      integer :: nom
      nom=num_omega1_adapt(i_nu)
      call zgemm('n','n',n,2*ndim3_tau,nom,(1.d0,0.d0),a_om,n,
     &           tau_from_omega_nu(1,1,1,i_nu),n_omega1_max,(0.d0,0.d0),
     &           a_tau,n)
      if(key==1) a_tau=conjg(a_tau)
      end



      subroutine from_omega_nu_to_tau_nu_aa(n,ind_nu,a_om,a_tau,key)
c   ------- Transform F(w;w-v) ---> A(t;v), B(t;v) --------
c     key=1   A
c     key=2   B
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,ind_nu,key
      complex*16, intent(in) :: a_om(n,n_omega1_max)
      complex*16, intent(out) :: a_tau(n,0:n_tau)
      integer :: nom,i_nu
	i_nu=me_t*ndim3_nu+ind_nu-1
      nom=num_omega1_adapt(i_nu)
      call zgemm('n','n',n,n_tau+1,nom,(1.d0,0.d0),a_om,n,
     &	       tau_nu_from_omega_nu(1,0,ind_nu),n_omega1_max,
     &           (0.d0,0.d0),a_tau,n)
      if(key==1) a_tau=conjg(a_tau)
      end
