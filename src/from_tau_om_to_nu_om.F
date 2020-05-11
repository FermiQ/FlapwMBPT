      subroutine from_tau_om_to_nu_om(n,ind_omega,nomnu,f_nu,f_tau)
c   ------- Transform AQ(t;v);BQ(t;v) ---> Q(v';v'-v) ---- N elements ----
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,ind_omega,nomnu
      real*8, intent(in) :: f_tau(n,0:n_tau/2,4)
      real*8, intent(out) :: f_nu(n,nomnu)
c ---------------------- COS(v')t contribution ---------------------------
      call dgemm('n','n',n,nomnu,n_tau/2+1,1.d0,f_tau(1,0,1),n,
     &           nu_om_from_tau_om(0,1,1,1,ind_omega),n_tau/2+1,0.d0,
     &           f_nu,n)
c ---------------------- SIN(v')t contribution ---------------------------
      call dgemm('n','n',n,nomnu,n_tau/2+1,1.d0,f_tau(1,0,2),n,
     &           nu_om_from_tau_om(0,1,2,1,ind_omega),n_tau/2+1,1.d0,
     &           f_nu,n)
c ---------------------- COS(v'-v)t contribution ---------------------------
      call dgemm('n','n',n,nomnu,n_tau/2+1,1.d0,f_tau(1,0,3),n,
     &           nu_om_from_tau_om(0,1,1,2,ind_omega),n_tau/2+1,1.d0,
     &           f_nu,n)
c ---------------------- SIN(v'-v)t contribution ---------------------------
      call dgemm('n','n',n,nomnu,n_tau/2+1,1.d0,f_tau(1,0,4),n,
     &           nu_om_from_tau_om(0,1,2,2,ind_omega),n_tau/2+1,1.d0,
     &           f_nu,n)
      end