      subroutine resp_uniform_prep_00_zz(p)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  real*8, intent(out) :: p
      integer :: k,ispin,i,j,ind_k,ind_tau,n,i_tau
      real*8 :: a,b
      complex*16, allocatable :: gt(:,:,:)
      p=0.d0
      do ispin=1,nspin
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          allocate(gt(n,n,2))
          b=0.d0
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            call ferm_unpack_tau(gt,g_full(1,1,1,ind_tau,ind_k,ispin),
     &	                         n,n,nbndf,3)
            a=0.d0
            do j=1,n
              do i=1,n
                a=a-gt(i,j,1)*gt(j,i,2)
              enddo
            enddo
            b=b+ab_nu_from_tau(i_tau,0,1)*a
          enddo
          p=p+wgt(k)*b
          deallocate(gt)
        enddo
      enddo
	  if(nproc_k/=1) call dgop(p,1,'  +',comm_k)
	  if(nproc_t/=1) call dgop(p,1,'  +',comm_t)
      p=2.d0*p/amega
      if(nspin==1.and.irel/=2) p=2.d0*p
      end
