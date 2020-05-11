      subroutine w_c_from_b
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ind_k,k,n,ind_nu,i_nu
      w_red_q=(0.d0,0.d0)
	  do ind_k=1,ndim3_k(me_k+1)
	    k=n3_mpi_k(me_k+1)+ind_k
        n=n_pbmt_red+nplwgw_red(k)
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          call unpack_hermit(w_red_q(1,1,ind_nu,k),
     &                       w_nu_q_b(1,1,ind_nu,ind_k),n,n_pbtot_red,
     &                       n_pbtot_red,(0.d0,0.d0),(1.d0,0.d0))
        enddo
      enddo
      if(nproc_k/=1) call dgop(w_red_q,2*n_pbtot_red**2*ndim3_nu*npnt,
     &                         '  +',comm_k)
      end
