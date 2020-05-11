      subroutine from_nu_omega_to_tau_omega_conjg(n,ind_omega,fnu,ftau)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,n
      complex*16, intent(in) :: fnu(n,n_nu1_max)
	complex*16, intent(out) :: ftau(n,0:n_tau,2)
      integer :: nn,n0,i_tau,i_omega
      real*8 :: phase
      complex*16 :: cc
	i_omega=me_t*ndim3_omega+ind_omega-1
c -------------------- A - component ( v <= 0 ) ------------------
      n0=nu_adapt_index_zero(i_omega)
      nn=n0	
	call zgemm('n','n',n,n_tau+1,nn,(1.d0,0.d0),fnu,n,
     &	         tau_from_nu_omega_1(1,0,ind_omega),n_nu1_max,
     &           (0.d0,0.d0),ftau(1,0,1),n)
c -------------------- B - component ( v >= 0 ) ------------------
      nn=num_nu1_adapt(i_omega)-n0+1
	call zgemm('n','n',n,n_tau+1,nn,(1.d0,0.d0),fnu(1,n0),n,
     &	         tau_from_nu_omega_1(n0,0,ind_omega),n_nu1_max,
     &           (0.d0,0.d0),ftau(1,0,2),n)
	do i_tau=0,n_tau
        phase=w_omega(i_omega)*tau_mesh(i_tau)
        cc=dcmplx(cos(phase),sin(phase))
        ftau(:,i_tau,2)=cc*ftau(:,i_tau,2)
      enddo
      end