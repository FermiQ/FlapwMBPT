      subroutine from_tau_nu_to_omega_nu_spl(n,ind_nu,nom,f_om,f_tau,
     &                                       key)
c   ------- Transform F(t;v) ---> F(-w+v;-w) : KEY = 1 -------
c   ------- Transform F(t;v) ---> F(w;w-v)  : KEY = 2 ------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,ind_nu,nom,key
      complex*16, intent(in) :: f_tau(n,0:n_tau/2,4)
      complex*16, intent(out) :: f_om(n,nom)
      complex*16, allocatable :: c(:,:),a(:,:)
      allocate(a(n,nom))
      allocate(c(n,nom))
      if(key==2) then   ! for w >= v/2
        call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,f_tau(1,0,1),2*n,
     &	         omega_nu_from_tau_nu_ab(0,1,1,ind_nu),n_tau/2+1,
     &             0.d0,a,2*n)
        call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,f_tau(1,0,2),2*n,
     &	         omega_nu_from_tau_nu_ab(0,1,2,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        a=a+(0.d0,1.d0)*c
        call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,f_tau(1,0,3),2*n,
     &             omega_nu_from_tau_nu_ab(0,1,3,ind_nu),n_tau/2+1,
     &             1.d0,a,2*n)
        call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,f_tau(1,0,4),2*n,
     &	         omega_nu_from_tau_nu_ab(0,1,4,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        a=a+(0.d0,1.d0)*c
      else if(key==1) then   ! for w <= v/2
        call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,f_tau(1,0,1),2*n,
     &	         omega_nu_from_tau_nu_ab(0,1,3,ind_nu),n_tau/2+1,
     &             0.d0,a,2*n)
        call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,f_tau(1,0,2),2*n,
     &          omega_nu_from_tau_nu_ab(0,1,4,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        a=a+(0.d0,-1.d0)*c
        call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,f_tau(1,0,3),2*n,
     &             omega_nu_from_tau_nu_ab(0,1,1,ind_nu),n_tau/2+1,
     &             1.d0,a,2*n)
        call dgemm('n','n',2*n,nom,n_tau/2+1,1.d0,f_tau(1,0,4),2*n,
     &	         omega_nu_from_tau_nu_ab(0,1,2,ind_nu),n_tau/2+1,
     &             0.d0,c,2*n)
        a=a+(0.d0,-1.d0)*c
      endif
      f_om=a
      deallocate(a,c)
      end