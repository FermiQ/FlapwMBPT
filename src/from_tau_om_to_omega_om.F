      subroutine from_tau_om_to_omega_om(n,ind_omega,nomom,f_om,f_tau,
     &                                   trn)
c   ------- Transform AQ(t;v);BQ(t;v) ---> Q(v';v'-v) ---- N elements ----
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,ind_omega,nomom
      real*8, intent(in) :: f_tau(n,0:n_tau/2,4),
     &                      trn(0:n_tau/2,n_omega1_max,4)
      real*8, intent(out) :: f_om(n,nomom)
c ---------------------- COS(w'+/-pi/b)t contribution --------------------
      call dgemm('n','n',n,nomom,n_tau/2+1,1.d0,f_tau(1,0,1),n,
     &           trn(0,1,1),n_tau/2+1,0.d0,f_om,n)
c ---------------------- SIN(w'+/-pi/b)t contribution --------------------
      call dgemm('n','n',n,nomom,n_tau/2+1,1.d0,f_tau(1,0,2),n,
     &           trn(0,1,2),n_tau/2+1,1.d0,f_om,n)
c ---------------------- COS(w-w')t contribution -----------------------
      call dgemm('n','n',n,nomom,n_tau/2+1,1.d0,f_tau(1,0,3),n,
     &           trn(0,1,3),n_tau/2+1,1.d0,f_om,n)
c ---------------------- SIN(w-w')t contribution -----------------------
      call dgemm('n','n',n,nomom,n_tau/2+1,1.d0,f_tau(1,0,4),n,
     &           trn(0,1,4),n_tau/2+1,1.d0,f_om,n)
      end