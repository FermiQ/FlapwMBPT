      subroutine set_g_full_0
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer :: ind_k,ispin,ln,k
      allocate(g_full_0(nbndf,nbndf,ndim3_k(me_k+1),nspin))
      g_full_0=(0.d0,0.d0)
      if(me_t==0) then
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            call ferm_unpack_tau(g_full_0(1,1,ind_k,ispin),
     &                           g_full(1,1,1,1,ind_k,ispin),
     &                           n_bnd(k,ispin),nbndf,nbndf,2)
            g_full_0(:,:,ind_k,ispin)=-g_full_0(:,:,ind_k,ispin)
          enddo   !! over k
        enddo  !! over ispin
      endif
      if(nproc_t/=1) then
        ln=nbndf**2*ndim3_k(me_k+1)*nspin
        if(ln/=0) call dgop(g_full_0,2*ln,'  +',comm_t)
      endif
      end
