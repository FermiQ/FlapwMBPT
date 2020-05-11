      subroutine lambda_bos_tau_matrix_r(lambda_store,lambda_matr,n0f,
     &                                   n0b,n1f,nf,jb)
c ---- We calculate dynamic tau-matrix for one fermionic index fixed ---- 
c --- TAU-matrix is:  LAM(t;t')=LAM(t;t')     -------------
c --- We don't need any interpolation here ------------------------
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: n0f,n0b,n1f,nf,jb
      real*8, intent(in) :: lambda_store(n0f*n0f,n0b,0:n_tau,2,
     &                                   ndim3_tau)
      real*8, intent(out) :: lambda_matr(n1f,n1f,0:n_tau,0:n_tau)
      integer :: i_tau,i1_tau,i,ind1_tau,it,j,ij
c ------------------------------------------------------------------
      lambda_matr=0.d0
      do ind1_tau=1,ndim3_tau
	  do it=1,2
	    i1_tau=me_t*ndim3_tau+ind1_tau-1
	    if(it==2) i1_tau=n_tau-i1_tau
          do i_tau=0,n_tau
            do j=1,nf
              do i=1,nf
                ij=n0f*(j-1)+i
	          lambda_matr(i,j,i_tau,i1_tau)=
     &	            lambda_store(ij,jb,i_tau,it,ind1_tau)
              enddo
            enddo
	    enddo  !! over i_tau
	  enddo
	enddo
	if(nproc_t/=1) call dgop(lambda_matr,n1f**2*(n_tau+1)**2,'  +',
     &	                         comm_t)
      end