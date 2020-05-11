      subroutine from_tau_nu_to_omega_nu(n,i_nu,nom,f_om,f_tau,key)
c   ------- Transform F(t;v) ---> F(-w+v;-w) : KEY = 1 -------
c   ------- Transform F(t;v) ---> F(w;w-v)  : KEY = 2 ------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,key,i_nu,nom
      complex*16, intent(in) :: f_tau(n,2,ndim3_tau,2)
      complex*16, intent(out) :: f_om(n,nom)
      integer :: ind_tau
      complex*16, allocatable :: c(:,:),a(:,:,:),b(:,:,:),d(:,:)
      allocate(c(n,nom))
      allocate(d(n,nom))
      allocate(a(n,ndim3_tau,2))
      allocate(b(n,ndim3_tau,2))
      do ind_tau=1,ndim3_tau
        a(:,ind_tau,1)=f_tau(:,1,ind_tau,1)-f_tau(:,2,ind_tau,1)
        a(:,ind_tau,2)=f_tau(:,1,ind_tau,2)-f_tau(:,2,ind_tau,2)
        b(:,ind_tau,1)=f_tau(:,1,ind_tau,1)+f_tau(:,2,ind_tau,1)
        b(:,ind_tau,2)=f_tau(:,1,ind_tau,2)+f_tau(:,2,ind_tau,2)
      enddo
      if(key==2) then   ! for w >= v/2
        call dgemm('n','n',2*n,nom,ndim3_tau,1.d0,a(1,1,1),2*n,
     &             omega_nu_from_tau_nu(1,1,1,i_nu),ndim3_tau,0.d0,d,
     &             2*n)
        call dgemm('n','n',2*n,nom,ndim3_tau,1.d0,b(1,1,1),2*n,
     &             omega_nu_from_tau_nu(1,1,2,i_nu),ndim3_tau,0.d0,c,
     &             2*n)
        d=d+(0.d0,1.d0)*c
        call dgemm('n','n',2*n,nom,ndim3_tau,1.d0,a(1,1,2),2*n,
     &             omega_nu_from_tau_nu(1,1,3,i_nu),ndim3_tau,1.d0,d,
     &             2*n)
        call dgemm('n','n',2*n,nom,ndim3_tau,1.d0,b(1,1,2),2*n,
     &             omega_nu_from_tau_nu(1,1,4,i_nu),ndim3_tau,0.d0,c,
     &             2*n)
        d=d+(0.d0,1.d0)*c
      else if(key==1) then   ! for w <= v/2
        call dgemm('n','n',2*n,nom,ndim3_tau,1.d0,a(1,1,1),2*n,
     &             omega_nu_from_tau_nu(1,1,3,i_nu),ndim3_tau,0.d0,d,
     &             2*n)
        call dgemm('n','n',2*n,nom,ndim3_tau,1.d0,b(1,1,1),2*n,
     &             omega_nu_from_tau_nu(1,1,4,i_nu),ndim3_tau,0.d0,c,
     &             2*n)
        d=d+(0.d0,-1.d0)*c
        call dgemm('n','n',2*n,nom,ndim3_tau,1.d0,a(1,1,2),2*n,
     &             omega_nu_from_tau_nu(1,1,1,i_nu),ndim3_tau,1.d0,d,
     &             2*n)
        call dgemm('n','n',2*n,nom,ndim3_tau,1.d0,b(1,1,2),2*n,
     &             omega_nu_from_tau_nu(1,1,2,i_nu),ndim3_tau,0.d0,c,
     &             2*n)
        d=d+(0.d0,-1.d0)*c
      endif
      f_om=d
      deallocate(c,a,b,d)
      if(nproc_t/=1) call dgop(f_om,2*n*nom,'  +',comm_t)
      end
