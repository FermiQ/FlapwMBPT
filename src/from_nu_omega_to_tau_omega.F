      subroutine from_nu_omega_to_tau_omega(i_omega,fnu,ftau,n)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: i_omega,n
      complex*16, intent(in) :: fnu(n,n_nu1_max)
	complex*16, intent(out) :: ftau(n,2,ndim3_tau,2)
      integer :: nn,n0,i_tau,it,ind_tau
      real*8 :: phase
      complex*16 :: cc
c -------------------- A - component ( v <= 0 ) ------------------
      n0=nu_adapt_index_zero(i_omega)
      nn=n0	
	call zgemm('n','n',n,2*ndim3_tau,nn,(1.d0,0.d0),fnu,n,
     &	         tau_from_nu_omega(1,1,1,i_omega),n_nu1_max,(0.d0,0.d0),
     &           ftau(1,1,1,1),n)
c -------------------- B - component ( v >= 0 ) ------------------
      nn=num_nu1_adapt(i_omega)-n0+1
	call zgemm('n','n',n,2*ndim3_tau,nn,(1.d0,0.d0),fnu(1,n0),n,
     &	         tau_from_nu_omega(n0,1,1,i_omega),n_nu1_max,
     &           (0.d0,0.d0),ftau(1,1,1,2),n)
	do ind_tau=1,ndim3_tau
	  do it=1,2
	    i_tau=me_t*ndim3_tau+ind_tau-1
	    if(it==2) i_tau=n_tau-i_tau
          phase=-w_omega(i_omega)*tau_mesh(i_tau)
          cc=dcmplx(cos(phase),sin(phase))
          ftau(:,it,ind_tau,2)=cc*ftau(:,it,ind_tau,2)
        enddo
      enddo
      end