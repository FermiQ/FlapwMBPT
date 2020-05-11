      subroutine a_cont_weights(achain,weigd,e,ispin,ind_k,n)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,ispin,ind_k
      real*8, intent(in) :: e
      complex*16, intent(in) :: achain(0:n_omega,n,n)
      real*8, intent(out) :: weigd(n)
      integer :: j,i,k,info
      complex*16 :: om0
      integer, allocatable :: ipi(:)
      complex*16, allocatable :: a(:,:),g_om(:,:)
      weigd=0.d0
      om0=dcmplx(e,e_small)
      allocate(ipi(n))
      allocate(a(n,n))
      allocate(g_om(n,n))
      k=n3_mpi_k(me_k+1)+ind_k
c ------ Chain approximation ----------------------------------
      call pade_fermion_all_out(n,n,achain,om0,a)
      g_om=(0.d0,0.d0)
      do j=1,n
        do i=1,n
          a(i,j)=-a(i,j)
          if(j==i) then
            a(i,j)=a(i,j)+om0+chem_pot-e_bnd(i,k,ispin)
            g_om(i,j)=(1.d0,0.d0)
          endif
        enddo
      enddo
      CALL ZGETRF(n,n,a,n,ipi,INFO)
      CALL ZGETRS('n',n,n,a,n,ipi,g_om,n,info)
      do i=1,n
        weigd(i)=-imag(g_om(i,i))/pi*wgt(k)
      enddo
      deallocate(a,g_om,ipi)
      weigd=weigd*2.d0/dfloat(nspin)/nrel
      end
