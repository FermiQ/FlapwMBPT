      subroutine from_omega_nu_to_tau_nu_d(n,i_nu,a_om,a_tau,key)
c   ------- Transform F(w;w-v) ---> F(t;v) --------
c     key = 1 : only for IT=1
c     key = 2 : only for IT=2
c     key = 3 : for both IT=1 and IT=2
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
	  if(key==1) then
	    call zgemm('n','n',n,ndim3_tau,nom,(1.d0,0.d0),a_om,n,
     &	           tau_from_omega_nu(1,1,1,i_nu),2*n_omega1_max,
     &             (0.d0,0.d0),a_tau(1,1,1),2*n)
	  else if(key==2) then
	    call zgemm('n','n',n,ndim3_tau,nom,(1.d0,0.d0),a_om,n,
     &	           tau_from_omega_nu(1,2,1,i_nu),2*n_omega1_max,
     &             (0.d0,0.d0),a_tau(1,2,1),2*n)
	  else if(key==3) then
	    call zgemm('n','n',n,2*ndim3_tau,nom,(1.d0,0.d0),a_om,n,
     &	           tau_from_omega_nu(1,1,1,i_nu),n_omega1_max,
     &             (0.d0,0.d0),a_tau(1,1,1),n)
      endif
      end
