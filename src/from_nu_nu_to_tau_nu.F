      subroutine from_nu_nu_to_tau_nu(n,ind_nu,f_nu,a_tau,b_tau)
c   ------- Transform K(v';v) ----> K(t;v) --------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n,ind_nu
	complex*16, intent(in) :: f_nu(n,n,n,n,n_nu2_max)
	complex*16, intent(out) :: a_tau(n,n,n,n,2,ndim3_tau),
     &	                         b_tau(n,n,n,n,2,ndim3_tau)
	integer :: i_nu,nom,j,i,n4,l,k
	i_nu=me_t*ndim3_nu+ind_nu-1
	nom=num_nu2_adapt(i_nu)
	n4=n**4
	call zgemm('n','n',n4,2*ndim3_tau,nom,(1.d0,0.d0),f_nu,n4,
     &	         tau_nu_from_nu_nu(1,1,1,ind_nu),n_nu2_max,(0.d0,0.d0),
     &           b_tau,n4)
      do l=1,n
        do k=1,n
          do j=1,n
            do i=1,n
              a_tau(i,j,k,l,:,:)=conjg(b_tau(j,i,l,k,:,:))
            enddo
          enddo
        enddo
      enddo
      end