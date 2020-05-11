      subroutine a_cont_coeff(achain,ispin,ind_k,n)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ispin,ind_k,n
      complex*16, intent(out) :: achain(0:n_omega,n,n)
      integer :: k,ind_omega,i_omega
      real*8 :: pib
      real*8, allocatable :: sigc(:,:,:,:)
      complex*16, allocatable :: sig_om(:,:,:)
      pib=pi/betta_t
      allocate(sigc(n,n,2,ndim3_omega))
      k=n3_mpi_k(me_k+1)+ind_k
      call tau_to_omega_baa(sigc,n,sig_c_tau(1,1,1,1,ind_k,ispin),nbndf,
     &                      n)
      allocate(sig_om(n,n,0:n_omega))
      sig_om=(0.d0,0.d0)
      do ind_omega=1,ndim3_omega
        i_omega=me_t*ndim3_omega+ind_omega-1
        call ferm_unpack_omega(sig_om(1,1,i_omega),
     &                         sigc(1,1,1,ind_omega),n,n,n)
      enddo
      deallocate(sigc)
      if(nproc_t/=1) call dgop(sig_om,2*n**2*(n_omega+1),'  +',comm_t)
      call pade_fermion_all(sig_om,n,n,achain)
      deallocate(sig_om)
      end
