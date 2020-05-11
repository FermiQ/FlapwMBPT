      subroutine sym_sig_k_red(s_bnd,s_rr,z_red,xx)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(in) :: s_rr(nrs_red,nrs_red,nqdiv_c,0:n_tau,
     &	                           nspin),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          xx(nr_full_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(out) :: s_bnd(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &	                             ndim_k_red(me_k+1),nspin)
      integer :: n,m,ispin,k0,nk,k,ind_tau,it,i_tau,ig,ig1,ind_k,ka0
      complex*16, allocatable :: s_in_mm(:,:),s_out_mm(:,:),s(:,:),
     &	                       tm(:,:),tt(:,:),s_in_mi(:,:),
     &                           s_out_mi(:,:),s_in_im(:,:),
     &                           s_out_im(:,:),s_in_ii(:,:),
     &                           s_out_ii(:,:)
      n=nfun_red
      m=nrs_red
      allocate(s(m,m))
      allocate(s_in_mm(n,n))
      allocate(s_out_mm(n,n))
      allocate(s_in_mi(n,m-n))
      allocate(s_out_mi(n,m-n))
      allocate(s_in_im(m-n,n))
      allocate(s_out_im(m-n,n))
      allocate(s_in_ii(m-n,m-n))
      allocate(s_out_ii(m-n,m-n))
      allocate(tm(m,nbndf_bnd))
      allocate(tt(m,nbndf_bnd))
      do ispin=1,nspin
        do ind_k=1,ndim_k_red(me_k+1)
          k0=n_mpi_k_red(me_k+1)+ind_k
          ka0=k_a_from_c(k0)
          nk=n_low_bnd(ka0,ispin)
          tm(1:n,1:nk)=z_red(:,1:nk,k0,ispin)
          tm(n+1:m,1:nk)=xx(:,1:nk,k0,ispin)
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=me_t*ndim3_tau+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
c ---------- Symmetrization -------------------------------------
              s=(0.d0,0.d0)
              do ig=1,ngroup
                ig1=u_inv(ig)  !! inverse op
                k=ip_k_c(k0,ig)  !! k=A^(-1)*k0
c ----------------- MT-MT part ----------------------------------
                s_in_mm=s_rr(1:n,1:n,k,i_tau,ispin)
                call sym_g_red_k(s_out_mm,s_in_mm,k0,ig)
                s(1:n,1:n)=s(1:n,1:n)+s_out_mm
c ----------------- MT-INT part ----------------------------------
                s_in_mi=s_rr(1:n,n+1:m,k,i_tau,ispin)
                call sym_g_red_k_mi(s_out_mi,s_in_mi,k0,ig)
                s(1:n,n+1:m)=s(1:n,n+1:m)+s_out_mi
c ----------------- INT-MT part ----------------------------------
                s_in_im=s_rr(n+1:m,1:n,k,i_tau,ispin)
                call sym_g_red_k_im(s_out_im,s_in_im,k0,ig)
                s(n+1:m,1:n)=s(n+1:m,1:n)+s_out_im
c ----------------- INT-INT part ----------------------------------
                s_in_ii=s_rr(n+1:m,n+1:m,k,i_tau,ispin)
                call sym_g_red_k_ii(s_out_ii,s_in_ii,k0,ig)
                s(n+1:m,n+1:m)=s(n+1:m,n+1:m)+s_out_ii
              enddo
              s=s/ngroup
c -------- Transform to BAND representation -----------------------
              call zgemm('n','n',m,nk,m,(1.d0,0.d0),s,m,tm,m,
     &	               (0.d0,0.d0),tt,m)
              call zgemm('c','n',nk,nk,m,(1.d0,0.d0),tm,m,tt,m,
     &	               (0.d0,0.d0),s_bnd(1,1,it,ind_tau,ind_k,ispin),
     &                   nbndf_bnd)
            enddo   !! over it
          enddo   !! over ind_tau
        enddo   !! over ind_k
      enddo    !! over ispin              
      deallocate(s,s_in_mm,s_out_mm,tm,tt,s_in_mi,s_out_mi,s_in_im,
     &           s_out_im,s_in_ii,s_out_ii)
      end
      
      subroutine sym_sig_k_red_pw(s_bnd,s_rr,z_red,xx,dd_fer)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(in) :: s_rr(nrs_red,nrs_red,nqdiv_c,0:n_tau,
     &	                           nspin),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          xx(nr_full_red,nbndf_bnd,nqdiv_c,nspin),
     &                          dd_fer(nr_full_red,nbndf_bnd,nqdiv_c,
     &                                 nspin)
      complex*16, intent(out) :: s_bnd(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &	                             ndim_k_red(me_k+1),nspin)
      integer :: n,m,ispin,k0,nk,k,ind_tau,it,i_tau,ig,ig1,ind_k,ka0
      complex*16 :: cc
      complex*16, allocatable :: s_in_mm(:,:),s_out_mm(:,:),s(:,:),
     &	                       tm(:,:),tt(:,:),s_in_mi(:,:),
     &                           s_out_mi(:,:),s_in_im(:,:),
     &                           s_out_im(:,:),s_in_ii(:,:),
     &                           s_out_ii(:,:)
      cc=amega/nr_full_red
      n=nfun_red
      m=nrs_red
      allocate(s(m,m))
      allocate(s_in_mm(n,n))
      allocate(s_out_mm(n,n))
      allocate(s_in_mi(n,m-n))
      allocate(s_out_mi(n,m-n))
      allocate(s_in_im(m-n,n))
      allocate(s_out_im(m-n,n))
      allocate(s_in_ii(m-n,m-n))
      allocate(s_out_ii(m-n,m-n))
      allocate(tm(m,nbndf_bnd))
      allocate(tt(nbndf_bnd,m))
      do ispin=1,nspin
        do ind_k=1,ndim_k_red(me_k+1)
          k0=n_mpi_k_red(me_k+1)+ind_k
          ka0=k_a_from_c(k0)
          nk=n_low_bnd(ka0,ispin)
          tm(1:n,1:nk)=z_red(:,1:nk,k0,ispin)
          tm(n+1:m,1:nk)=xx(:,1:nk,k0,ispin)
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=me_t*ndim3_tau+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
c ---------- Symmetrization -------------------------------------
              s=(0.d0,0.d0)
              do ig=1,ngroup
                ig1=u_inv(ig)  !! inverse op
                k=ip_k_c(k0,ig)  !! k=A^(-1)*k0
c ----------------- MT-MT part ----------------------------------
                s_in_mm=s_rr(1:n,1:n,k,i_tau,ispin)
                call sym_g_red_k(s_out_mm,s_in_mm,k0,ig)
                s(1:n,1:n)=s(1:n,1:n)+s_out_mm
c ----------------- MT-INT part ----------------------------------
                s_in_mi=s_rr(1:n,n+1:m,k,i_tau,ispin)
                call sym_g_red_k_mi(s_out_mi,s_in_mi,k0,ig)
                s(1:n,n+1:m)=s(1:n,n+1:m)+s_out_mi
c ----------------- INT-MT part ----------------------------------
                s_in_im=s_rr(n+1:m,1:n,k,i_tau,ispin)
                call sym_g_red_k_im(s_out_im,s_in_im,k0,ig)
                s(n+1:m,1:n)=s(n+1:m,1:n)+s_out_im
c ----------------- INT-INT part ----------------------------------
                s_in_ii=s_rr(n+1:m,n+1:m,k,i_tau,ispin)
                call sym_g_red_k_ii(s_out_ii,s_in_ii,k0,ig)
                s(n+1:m,n+1:m)=s(n+1:m,n+1:m)+s_out_ii
              enddo
              s=s/ngroup
c -------- Transform to BAND representation -----------------------
              call zgemm('c','n',nk,m,m,(1.d0,0.d0),tm,m,s,m,
     &	               (0.d0,0.d0),tt,nbndf_bnd)
              call zgemm('n','n',nk,nk,nfun_red,(1.d0,0.d0),tt,
     &                   nbndf_bnd,z_red(1,1,k0,ispin),nfun_red,
     &                   (0.d0,0.d0),s_bnd(1,1,it,ind_tau,ind_k,ispin),
     &                   nbndf_bnd)
              call zgemm('n','n',nk,nk,nr_full_red,cc,tt(1,nfun_red+1),
     &                   nbndf_bnd,dd_fer(1,1,k0,ispin),nr_full_red,
     &                   (1.d0,0.d0),s_bnd(1,1,it,ind_tau,ind_k,ispin),
     &                   nbndf_bnd)
            enddo   !! over it
          enddo   !! over ind_tau
        enddo   !! over ind_k
      enddo    !! over ispin              
      deallocate(s,s_in_mm,s_out_mm,tm,tt,s_in_mi,s_out_mi,s_in_im,
     &           s_out_im,s_in_ii,s_out_ii)
      end