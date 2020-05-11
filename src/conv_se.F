      subroutine conv_se(conv_vh,conv_sx,conv_sc)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(out) :: conv_vh,conv_sx,conv_sc
      integer :: ispin,ind_k,k,i,j,ind_tau,it,n,nn,ind_ij,ij
      real*8 :: x
      conv_vh=0.d0
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          nn=n*(n+1)/2
          do ind_ij=1,nd_tb_ij(me_t_b+1)
            ij=n0_tb_ij(me_t_b+1)+ind_ij
            if(ij>nn) cycle
            x=abs(v_h_bb(ind_ij,ind_k,ispin)
     &           -v_h_bb_old(ind_ij,ind_k,ispin))
            conv_vh=max(conv_vh,x)
          enddo
        enddo
      enddo
      if(nproc_k/=1) call dgop(conv_vh,1,'max',comm_k)
      if(nproc_t*nproc_k/=1) call dgop(conv_vh,1,'max',comm_t_b)
c ------------------------------------------------------------------
      conv_sx=0.d0
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          nn=n*(n+1)/2
          do ind_ij=1,nd_tb_ij(me_t_b+1)
            ij=n0_tb_ij(me_t_b+1)+ind_ij
            if(ij>nn) cycle
            conv_sx=max(conv_sx,abs(sigx_solid_k(ind_ij,ind_k,ispin)
     &                          -sigx_solid_k_old(ind_ij,ind_k,ispin)))
          enddo
        enddo
      enddo
      if(nproc_k/=1) call dgop(conv_sx,1,'max',comm_k)
      if(nproc_t*nproc_k/=1) call dgop(conv_sx,1,'max',comm_t_b)
      if(ubi==' hf') return
c ------------------------------------------------------------------
      conv_sc=0.d0
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          do ind_tau=1,ndim3_tau
            do it=1,2
              do j=1,n
                do i=1,n
                  conv_sc=max(conv_sc,
     &                        abs(sig_c_tau(i,j,it,ind_tau,ind_k,ispin)
     &                     -sig_c_tau_old(i,j,it,ind_tau,ind_k,ispin)))
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      if(nproc_k/=1) call dgop(conv_sc,1,'max',comm_k)
      if(nproc_t/=1) call dgop(conv_sc,1,'max',comm_t)
      end
