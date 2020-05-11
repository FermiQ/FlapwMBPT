      subroutine sym_sig_k_red_a(s_bnd,s_rr,z_red,a_red)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	complex*16, intent(in) :: s_rr(nra_red,nra_red,nqdiv_c,0:n_tau,
     &	                           nspin),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,
     &                                nspin)
	complex*16, intent(out) :: s_bnd(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &	                             ndim_k_red(me_k+1),nspin)
	integer :: n,m,ispin,k0,nk,k,ind_tau,it,i_tau,ig,ind_k,ka0,iga,m0
	complex*16, allocatable :: s_in_mm(:,:),s_out_mm(:,:),s(:,:),
     &	                       tm(:,:),tt(:,:),s_in_mi(:,:),
     &                           s_out_mi(:,:),s_in_im(:,:),
     &                           s_out_im(:,:),s_in_ii(:,:),
     &                           s_out_ii(:,:),st(:,:,:)
	n=nfun_red
	m=nra_red
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
          m0=n+nbask_red(ka0)
	    tm(1:n,1:nk)=z_red(:,1:nk,k0,ispin)
	    tm(n+1:m0,1:nk)=a_red(1:m0-n,1:nk,k0,ispin)
          allocate(st(nk,nk,0:n_tau))
          st=(0.d0,0.d0)
	    do ind_tau=1,ndim3_tau
	      do it=1,2
	        i_tau=me_t*ndim3_tau+ind_tau-1
	        if(it==2) i_tau=n_tau-i_tau
c ---------- Symmetrization -------------------------------------
              s=(0.d0,0.d0)
              do ig=1,ngroup
                k=k_sym_0_c(k0,ig)  !! k=A^(-1)*k0
                iga=g_sym_0_c(k0,ig)
c ----------------- MT-MT part ----------------------------------
                s_in_mm=s_rr(1:n,1:n,k,i_tau,ispin)
                call sym_g_red_k(s_out_mm,s_in_mm,k0,ig)
                s(1:n,1:n)=s(1:n,1:n)+s_out_mm
c ----------------- MT-INT part ----------------------------------
                s_in_mi=s_rr(1:n,n+1:m,k,i_tau,ispin)
                call sym_sig_red_k_mg(s_out_mi,s_in_mi,k0,k,ig,iga)
                s(1:n,n+1:m)=s(1:n,n+1:m)+s_out_mi
c ----------------- INT-MT part ----------------------------------
                s_in_im=s_rr(n+1:m,1:n,k,i_tau,ispin)
                call sym_sig_red_k_gm(s_out_im,s_in_im,k0,k,ig,iga)
                s(n+1:m,1:n)=s(n+1:m,1:n)+s_out_im
c ----------------- INT-INT part ----------------------------------
                s_in_ii=s_rr(n+1:m,n+1:m,k,i_tau,ispin)
                call sym_sig_red_k_gg(s_out_ii,s_in_ii,k0,k,ig,iga)
                s(n+1:m,n+1:m)=s(n+1:m,n+1:m)+s_out_ii
              enddo
              s=s/ngroup
c -------- Transform to BAND representation -----------------------
	        call zgemm('n','n',m0,nk,m0,(1.d0,0.d0),s,m,tm,m,
     &	               (0.d0,0.d0),tt,m)
	        call zgemm('c','n',nk,nk,m0,(1.d0,0.d0),tm,m,tt,m,
     &	               (0.d0,0.d0),st(1,1,i_tau),nk)
            enddo   !! over it
          enddo   !! over ind_tau
	    if(nproc_t/=1) call dgop(st,2*nk*nk*(n_tau+1),'  +',
     &                               comm_t)
	    do ind_tau=1,ndim3_tau
	      do it=1,2
	        i_tau=me_t*ndim3_tau+ind_tau-1
	        if(it==2) i_tau=n_tau-i_tau
              s_bnd(1:nk,1:nk,it,ind_tau,ind_k,ispin)=st(:,:,i_tau)
            enddo
          enddo
          deallocate(st)
        enddo   !! over ind_k
      enddo    !! over ispin              
      deallocate(s,s_in_mm,s_out_mm,tm,tt,s_in_mi,s_out_mi,s_in_im,
     &           s_out_im,s_in_ii,s_out_ii)
      end
      
      
      subroutine sym_sig_k_red_0(s_bnd,s_rr,z_red)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	complex*16, intent(in) :: s_rr(nfun_red,nfun_red,nqdiv_c,0:n_tau,
     &	                           nspin),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
	complex*16, intent(out) :: s_bnd(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &	                             ndim_k_red(me_k+1),nspin)
	integer :: n,ispin,k0,nk,k,ind_tau,it,i_tau,ig,ind_k,ka0,iga
	complex*16, allocatable :: s_in_mm(:,:),s_out_mm(:,:),s(:,:),
     &	                       tm(:,:),tt(:,:),st(:,:,:)
	n=nfun_red
	allocate(s(n,n))
	allocate(s_in_mm(n,n))
	allocate(s_out_mm(n,n))
	allocate(tm(n,nbndf_bnd))
	allocate(tt(n,nbndf_bnd))
	do ispin=1,nspin
        do ind_k=1,ndim_k_red(me_k+1)
          k0=n_mpi_k_red(me_k+1)+ind_k
          ka0=k_a_from_c(k0)
	    nk=n_low_bnd(ka0,ispin)
	    tm(1:n,1:nk)=z_red(:,1:nk,k0,ispin)
          allocate(st(nk,nk,0:n_tau))
          st=(0.d0,0.d0)
	    do ind_tau=1,ndim3_tau
	      do it=1,2
	        i_tau=me_t*ndim3_tau+ind_tau-1
	        if(it==2) i_tau=n_tau-i_tau
c ---------- Symmetrization -------------------------------------
              s=(0.d0,0.d0)
              do ig=1,ngroup
                k=k_sym_0_c(k0,ig)  !! k=A^(-1)*k0
                iga=g_sym_0_c(k0,ig)
c ----------------- MT-MT part ----------------------------------
                s_in_mm=s_rr(1:n,1:n,k,i_tau,ispin)
                call sym_g_red_k(s_out_mm,s_in_mm,k0,ig)
                s(1:n,1:n)=s(1:n,1:n)+s_out_mm
              enddo
              s=s/ngroup
c -------- Transform to BAND representation -----------------------
	        call zgemm('n','n',n,nk,n,(1.d0,0.d0),s,n,tm,n,
     &	               (0.d0,0.d0),tt,n)
	        call zgemm('c','n',nk,nk,n,(1.d0,0.d0),tm,n,tt,n,
     &	               (0.d0,0.d0),st(1,1,i_tau),nk)
            enddo   !! over it
          enddo   !! over ind_tau
	    if(nproc_t/=1) call dgop(st,2*nk*nk*(n_tau+1),'  +',
     &                               comm_t)
	    do ind_tau=1,ndim3_tau
	      do it=1,2
	        i_tau=me_t*ndim3_tau+ind_tau-1
	        if(it==2) i_tau=n_tau-i_tau
              s_bnd(1:nk,1:nk,it,ind_tau,ind_k,ispin)=st(:,:,i_tau)
            enddo
          enddo
          deallocate(st)
        enddo   !! over ind_k
      enddo    !! over ispin              
      deallocate(s,s_in_mm,s_out_mm,tm,tt)
      end