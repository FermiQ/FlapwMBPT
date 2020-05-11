      subroutine set_g00
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      logical :: correlation
      integer :: k,ispin,i,ind_k
      correlation=.false.
      if(ubi==' qp'.or.ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') 
     &                                           correlation=.true.
      if(ubi==' qp') then
        if(ubi_0/=' qp') then
          e_qp=e_bnd
          z_qp=z_bnd
          q_qp=(0.d0,0.d0)
          do ispin=1,nspin
            do ind_k=1,ndim3_k(me_k+1)
              k=n3_mpi_k(me_k+1)+ind_k
              do i=1,n_bnd(k,ispin)
                q_qp(i,i,ind_k,ispin)=(1.d0,0.d0)
              enddo
            enddo
          enddo
        endif
      endif
      if(ubi==' gw'.and.ubi_0/=' gw') then
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me_k+1)
            call set_g_tau_0(ispin,ind_k,0,
     &	                     g_full(1,1,1,1,ind_k,ispin))
          enddo
        enddo
      endif
      if(ubi==' qp'.and.ubi_0/=' qp'.and.ubi_0/=' gw') then
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me_k+1)
            call set_g_tau_0(ispin,ind_k,0,
     &	                     g_full(1,1,1,1,ind_k,ispin))
          enddo
        enddo
      endif
      if(ubi=='psi'.and.ubi_0/='psi'.and.ubi_0/=' gw'.
     &                               and.ubi_0/=' qp') then
        if(.not.allocated(g_full)) 
     &	  allocate(g_full(nbndf,nbndf,2,ndim3_tau,ndim3_k(me_k+1),
     &                      nspin))
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me_k+1)
            call set_g_tau_0(ispin,ind_k,0,
     &	                   g_full(1,1,1,1,ind_k,ispin))
          enddo
        enddo
      endif
      if(ubi=='bsp'.and.ubi_0/='bsp'.and.ubi_0/=' gw'.
     &                and.ubi_0/=' qp'.and.ubi_0/='psi') then
        if(.not.allocated(g_full)) 
     &	  allocate(g_full(nbndf,nbndf,2,ndim3_tau,ndim3_k(me_k+1),
     &                      nspin))
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me_k+1)
            call set_g_tau_0(ispin,ind_k,0,
     &	                   g_full(1,1,1,1,ind_k,ispin))
          enddo
        enddo
      endif
      if(iter/=1) then
c ------------- Forming G_FULL(0-) --> G_FULL_0 --------------------
        call green_0(correlation)
c ------------- Forming Local Green's Function ---------------------
        call g_mt_from_g_full(correlation,g_loc_0)
      endif
      end
