      subroutine sigc_r_tau
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: ispin,ind_k,k,ib,ind_tau,it,nbnd,jb,i_tau
      real*8 :: dum,w
      complex*16 :: dumc
      real*8, allocatable :: gr1(:,:,:),tmp2(:,:,:)
      complex*16, allocatable :: s_c1(:,:,:,:),gc1(:,:,:,:),
     &                           sig_tmp(:,:,:,:,:),tmp1(:,:,:)
c -------- Precalculation of the matrix S*A ------------------------
      allocate(sa(nplw_gw,nbndf,ndim3_k(me_k+1),nrel,nspin))
      call sa_prep(sa)
c ------------------------------------------------------------------
      allocate(sig_tmp(nbndf,nbndf,2,ndim3_k(me_k+1),nspin))
      do ind_tau=1,ndim3_tau
        sig_tmp=(0.d0,0.d0)
c
c	---------- MT_MT Part ------------------------------------
c
        call sig_c_mm(ind_tau)
        if(ind_tau==1)
     &          call timel('** SIG_C_MM finished ***************')
c
c	---------- INT_MT Part ------------------------------------
c		    
        allocate(w_rs_mi(nr_red,n_pbmt,ndim3_kk(me_k+1)))
        call w_im_r_from_k(ind_tau)
        if(ind_tau==1)
     &      call timel('** W_IM_R_FROM_K_B finished ********')
        allocate(s_c1(nr_red,nfun,nrel,ndim3_kk(me_k+1)))
        do ispin=1,nspin
          do it=1,2
            if(irel<=1) then
              allocate(gr1(nr_red,nfun,ndim3_kk(me_k+1)))
              call g_rs_mi_mem(gr1,dumc,it,ispin,ind_tau,0)
              if(ind_tau==1.and.it==1)
     &            call timel('** G_RS_FROM_KS_MI finished ********')
              call sigc_im_r(ispin,s_c1,gr1,dumc)
              deallocate(gr1)
            else if(irel==2) then
              allocate(gc1(nr_red,nfun,ndim3_kk(me_k+1),nrel))
              call g_rs_mi_mem(dum,gc1,it,1,ind_tau,0)
              if(ind_tau==1.and.it==1)
     &            call timel('** G_RS_FROM_KS_MI finished ********')
              call sigc_im_r(ispin,s_c1,dum,gc1)
              deallocate(gc1)
            endif
            if(ind_tau==1.and.it==1)
     &          call timel('** SIGC_IM_R finished **************')
            call sigc_im_r_to_k(ispin,it,s_c1,sig_tmp(1,1,1,1,ispin))
            if(ind_tau==1.and.it==1)
     &          call timel('** SIGC_IM_R_TO_K finished *********')
          enddo
        enddo  !! over ispin
        deallocate(w_rs_mi,s_c1)
c
c	---------- INT_INT Part ------------------------------------
c		  
        call sigc_ii_mem(ind_tau,sig_tmp)
        allocate(tmp2(nbndf,nbndf,2))
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            nbnd=n_bnd(k,ispin)
            tmp2=0.d0
            call ferm_pack_tau(sig_tmp(1,1,1,ind_k,ispin),tmp2,nbnd,
     &                         nbndf,nbndf)
            do jb=1,nbnd
              do ib=1,nbnd
                do it=1,2
                  sig_c_tau(ib,jb,it,ind_tau,ind_k,ispin)=
     &              sig_c_tau(ib,jb,it,ind_tau,ind_k,ispin)
     &             +tmp2(ib,jb,it)
                enddo
              enddo
            enddo
          enddo   !! over ind_k
        enddo
        deallocate(tmp2)
      enddo   !! over ind_tau
      deallocate(sig_tmp,sa)
      allocate(tmp1(nbndf,nbndf,2))
      w=-wgt(1)/amega*q2aver
      do ind_tau=1,ndim3_tau
        i_tau=me_t*ndim3_tau+ind_tau-1
        do ispin=1,nspin
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            nbnd=n_bnd(k,ispin)
            call ferm_unpack_tau(tmp1,g_full(1,1,1,ind_tau,ind_k,ispin),
     &                           nbnd,nbndf,nbndf,3)
            do jb=1,nbnd
              do ib=1,nbnd
                do it=1,2
                  tmp1(ib,jb,it)=w*ws_head(i_tau)*tmp1(ib,jb,it)
                enddo
              enddo
            enddo
            call ferm_pack_tau(tmp1,
     &                         sig_c_tau(1,1,1,ind_tau,ind_k,ispin),
     &                         nbnd,nbndf,nbndf)
          enddo   !! over ind_k
        enddo
      enddo   !! over ind_tau
      deallocate(tmp1)
      end
