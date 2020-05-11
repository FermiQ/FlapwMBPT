      subroutine bands
****************************************************************
*   Control module for package of programs solved              *
*         eigenvalue problem of the APW method                 *
****************************************************************
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer :: k,ind_k,ispin
      real*8 :: err_val,err_slo
* ****** compute eigenvalues and eigenvectors ************
      err_val=0.d0
      err_slo=0.d0
      if(l_grad_k) de_bnd=0.d0
      do ispin=1,nspin      !!! do over spins
        if(nspin.eq.2.and.maswrk)
     &    write(iun,*)'&&&&&&&&&   ISPIN   = ',ispin,'  &&&&&&&&&&&&'
        n_bnd(:,ispin)=0
        ev_bnd(:,:,:,ispin)=(0.d0,0.d0)
        z_bnd(:,:,:,ispin)=(0.d0,0.d0)
        e_bnd(:,:,ispin)=0.d0
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          call lda_k_point(ind_k,ispin,pnt(1,k),err_val,err_slo,
     &	                   e_bnd(1,k,ispin),n_bnd(k,ispin),
     &                     z_bnd(1,1,ind_k,ispin),
     &                     ev_bnd(1,1,ind_k,ispin),indgb(1,k),nbask(k),
     &                     0)
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          if(k.eq.1) call timel('****** LDA_K_POINT finished ********')
        enddo                !!! over k_ind
        if(nproc_k/=1) then
          call igop(n_bnd(1,ispin),npnt,'  +',comm_k)
          call dgop(e_bnd(1,1,ispin),nbndf*npnt,'  +',comm_k)
        endif
c -------------------------------------------------------------------
      enddo   !! over ispin
      if(nproc_k/=1) then
        call dgop(err_val,1,'max',comm_k)
        call dgop(err_slo,1,'max',comm_k)
      endif                !!! over k_ind
      if(l_grad_k.and.nproc_k/=1) call dgop(de_bnd,3*nbndf*npnt*nspin,
     &                                      '  +',comm_k)
      if(maswrk) write(iun,'(a39,2e12.5)')
     &	' BANDS: Max mismatch in value, slope : ',err_val,err_slo
c ------------------------------------------------------------------
      call bnd_k_blocks
      call bnd_curv
      call timel('******* bands finished *************')
      end
