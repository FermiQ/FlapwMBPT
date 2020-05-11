      subroutine vertex_wk_uni_test(fif_full,k_pw,g_tau,z_full,a_full,
     &                              w_tau,lambda_dyn,lambda_stat,v2_q)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: fif_full(nrel,maxel**2,n_pbmtm_red,nsort,
     &                               nspin),
     &                      k_pw(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                           npnt_c,nspin),
     &                      g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c,nspin),
     &                      w_tau(n_pbtot_red,n_pbtot_red,ndim3_tau,
     &                            nqdiv_c)
      complex*16, intent(in) :: z_full(nfun,nbndf_bnd,nqdiv_c,nspin),
     &                          a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin),
     &                          v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c)
      real*8, intent(out) :: lambda_dyn(nbndf_bnd,nbndf_bnd,2,
     &                                  ndim3_omega,npnt_c,nspin)
      complex*16, intent(out) :: lambda_stat(nbndf_bnd,nbndf_bnd,npnt_c,
     &                                       nspin)
      integer :: k,ispin,nn,n,i_tau,nk,i,ka,kq,kq0,kqa,nkq,iqa,iq0,
     &           ind_tau,iq,it,j,kk,l,np,ind_k,ind_q
      real*8 :: v(3),gtild(3,2),rc
      integer, allocatable :: nd(:),n0(:)
      real*8, allocatable :: vrt0(:,:,:,:)
      complex*16, allocatable :: tmp(:,:),ps3(:,:,:),a(:,:,:),
     &                           kt0(:,:),b(:,:,:,:),c(:,:,:,:),
     &                           wt0(:,:),vrt(:,:,:,:)
      n=nbndf_bnd
      nn=n*n
      allocate(ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red))
      allocate(kt0(n,n))
      allocate(tmp(n,n))
      allocate(nd(nproc_b))
      allocate(n0(nproc_b))
      call size_shift_par(nqdiv_c,nproc_b,nd,n0)
      lambda_stat=(0.d0,0.d0)
      lambda_dyn=0.d0
      do ispin=1,nspin
        do ind_k=1,ndim_k_red(me_k+1)
          k=n_mpi_k_red(me_k+1)+ind_k
          ka=k_a_from_c(k)
          nk=n_low_bnd(ka,ispin)
          allocate(vrt(nk,nk,2,ndim3_tau))
          vrt=(0.d0,0.d0)
          do ind_q=1,nd(me_b+1)
            iq=n0(me_b+1)+ind_q
            iq0=i_kref_c(iq)
            v=pnt_c(:,k)-pnt_c(:,iq)
            iqa=k_a_from_c(iq0)
            np=n_pbmt_red+nplwgw_red(iqa)
            call zone1_number(v,rb0_c,ndiv_c,kq)
            kq=index_k1_c(kq)    !! for K-Q
            kq0=i_kref_c(kq)
            kqa=k_a_from_c(kq0)
            nkq=n_low_bnd(kqa,ispin)
            gtild(:,1)=0.d0
            gtild(:,2)=v-pnt_c(:,kq)
            call ppm_factors_full(k,kq,iq,fif_full(1,1,1,1,ispin),
     &                            z_full(1,1,1,ispin),
     &                            a_full(1,1,1,ispin),ps3,gtild,
     &                            n_low_bnd(1,ispin))
            allocate(a(n,n,n_pbtot_red))
            allocate(b(n,n,n,n))
            allocate(c(n,n,n,n))
c --------------- Static part --------------------------------------
            call zgemm('n','n',nn,np,np,(1.d0,0.d0),ps3,nn,
     &                 v2_q(1,1,iq),n_pbtot_red,(0.d0,0.d0),a,nn)
            if(me_t==0) then
              call ferm_unpack_tau(kt0,k_pw(1,1,1,1,kq0,ispin),nkq,n,
     &                             nbndf_bnd,1)
            endif
            if(nproc_t/=1) call brdcst(msgdbl,kt0,16*n*n,0,comm_t)
            call zgemm('n','c',nn,nn,np,(1.d0,0.d0),a,nn,ps3,nn,
     &                 (0.d0,0.d0),b,nn)
            do l=1,n
              do kk=1,n
                do j=1,n
                  do i=1,n
                    c(i,j,kk,l)=b(i,kk,j,l)
                  enddo
                enddo
              enddo
            enddo
            call zgemm('n','n',nn,1,nn,(1.d0,0.d0),c,nn,kt0,nn,
     &                 (0.d0,0.d0),tmp,nn)
            do j=1,nk
              do i=1,nk
                lambda_stat(i,j,k,ispin)=lambda_stat(i,j,k,ispin)
     &                                  +tmp(i,j)/nqdiv_c
              enddo
            enddo
c ---------------- Dynamic part ------------------------------------
            allocate(wt0(np,np))
            do ind_tau=1,ndim3_tau
              call boson_unpack_tau(wt0,np,w_tau(1,1,ind_tau,iq),
     &                              n_pbtot_red,np)
              call zgemm('n','n',nn,np,np,(1.d0,0.d0),ps3,nn,wt0,np,
     &                   (0.d0,0.d0),a,nn)
              call zgemm('n','c',nn,nn,np,(1.d0,0.d0),a,nn,ps3,nn,
     &                   (0.d0,0.d0),b,nn)
              do l=1,n
                do kk=1,n
                  do j=1,n
                    do i=1,n
                      c(i,j,kk,l)=b(i,kk,j,l)
                    enddo
                  enddo
                enddo
              enddo
              do it=1,2
                call ferm_unpack_tau(kt0,k_pw(1,1,1,ind_tau,kq0,ispin),
     &                               nkq,nbndf_bnd,nbndf_bnd,it)
                call zgemm('n','n',nn,1,nn,(1.d0,0.d0),c,nn,kt0,nn,
     &                     (0.d0,0.d0),tmp,nn)
                do j=1,nk
                  do i=1,nk
                    vrt(i,j,it,ind_tau)=vrt(i,j,it,ind_tau)
     &                                 +tmp(i,j)/nqdiv_c
                  enddo
                enddo
              enddo
            enddo   !! over ind_tau
            deallocate(a,b,c,wt0)
          enddo   !! over iq
          if(nproc_b/=1) call dgop(vrt,4*nk**2*ndim3_tau,'  +',comm_b)
c ---------- Static Singular part ----------------------------------
          rc=8.d0*pi*q2aver_c/amega*wgt_c(1)
          if(me_t==0) then
            call ferm_unpack_tau(kt0,k_pw(1,1,1,1,k,ispin),nk,n,
     &                           nbndf_bnd,1)
          endif
          if(nproc_t/=1) call brdcst(msgdbl,kt0,16*n*n,0,comm_t)
          do j=1,nk
            do i=1,nk
              lambda_stat(i,j,k,ispin)=lambda_stat(i,j,k,ispin)
     &                                +rc*kt0(i,j)
            enddo
          enddo
c ---------------- Dynamic Singular part ---------------------------
          allocate(vrt0(nk,nk,2,ndim3_tau))
          vrt0=0.d0
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            rc=ws_head(i_tau)*q2aver_c/amega*wgt_c(1)
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              call ferm_unpack_tau(kt0,k_pw(1,1,1,ind_tau,k,ispin),nk,n,
     &                             nbndf_bnd,it)
              do j=1,nk
                do i=1,nk
                  vrt(i,j,it,ind_tau)=vrt(i,j,it,ind_tau)+rc*kt0(i,j)
                enddo
              enddo
            enddo
            call ferm_pack_tau(vrt(1,1,1,ind_tau),vrt0(1,1,1,ind_tau),
     &                         nk,nk,nk)
          enddo
          deallocate(vrt)
          call tau_to_omega_baa(lambda_dyn(1,1,1,1,k,ispin),nbndf_bnd,
     &                          vrt0,nk,nk)
          deallocate(vrt0)
        enddo   !! over k
      enddo   !! over ispin
      deallocate(ps3,kt0,tmp,nd,n0)
      if(nproc_k/=1) then
        call dgop(lambda_stat,2*nbndf_bnd**2*npnt_c*nspin,'  +',comm_k)
        call dgop(lambda_dyn,2*nbndf_bnd**2*ndim3_omega*npnt_c*nspin,
     &            '  +',comm_k)
      endif
      if(nproc_b/=1) then
        call dgop(lambda_stat,2*nbndf_bnd**2*npnt_c*nspin,'  +',comm_b)
      endif
      end
