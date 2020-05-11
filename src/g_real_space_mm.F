      subroutine g_real_space_mm(g_mm,ispin)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      real*8, intent(out) :: g_mm(nfun_red,nfun_red,nqdiv_c,2,ndim3_tau)
      integer :: k,ind_tau,it,nb0,i_tau,nb,k0,k1,k11,ir,kf,ind_k,kc,n,
     &           ibnd,ib,ig,jbnd,jb
      real*8, allocatable :: g_tau(:,:,:,:,:)
      complex*16, allocatable :: g(:,:),tm(:,:,:),tmp(:,:),tmp1(:,:),
     &                           tt(:,:)
	  nb0=nbndf_bnd
	  g_mm=0.d0
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c))
      g_tau=0.d0
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        kc=k_c_from_a(k)
        if(kc==0) cycle
        n=n_low_bnd(k,ispin)
        do ind_tau=1,ndim3_tau
          do it=1,2
	        do jbnd=1,n
	          jb=ind_bands_bnd(jbnd,k,ispin)
	          do ibnd=1,n
	            ib=ind_bands_bnd(ibnd,k,ispin)
                g_tau(ibnd,jbnd,it,ind_tau,kc)=
     &	            g_full(ib,jb,it,ind_tau,ind_k,ispin)
              enddo
            enddo
          enddo
        enddo   !! over ind_tau
	  enddo  ! over ind_k
	  if(nproc_k/=1) then
	    call dgop(g_tau,2*nbndf_bnd**2*npnt_c*ndim3_tau,'  +',
     &            comm_k)
      endif
	  allocate(g(nb0,nb0))
	  allocate(tmp(nfun_red,nb0))
	  allocate(tmp1(nfun_red,nfun_red))
	  allocate(tm(nqdiv_c,nfun_red,nfun_red))
      allocate(tt(nfun_red,nfun_red))
	  do ind_tau=1,ndim3_tau
        i_tau=ndim3_tau*me_t+ind_tau-1
        do it=1,2
          if(it==2) i_tau=n_tau-i_tau
          tm=(0.d0,0.d0)
	      do k=1,npnt_c
            kf=k_a_from_c(k)
	        nb=n_low_bnd(kf,ispin)
c -------------------- G on K_TAU ------------------------------------
            call ferm_unpack_tau(g,g_tau(1,1,1,ind_tau,k),nb,nb0,nb0,it)
	        call zgemm('n','n',nfun_red,nb,nb,(1.d0,0.d0),
     &                 z_bnd_red(1,1,kf,ispin),nfun_red,g,nb0,
     &                 (0.d0,0.d0),tmp,nfun_red)
            call zgemm('n','c',nfun_red,nfun_red,nb,(1.d0,0.d0),tmp,
     &	               nfun_red,z_bnd_red(1,1,kf,ispin),nfun_red,
     &                 (0.d0,0.d0),tmp1,nfun_red)
            do k0=1,k_star_c(k)
	          k1=k_list_c(k0,k)
              ig=k_group_c(k1)
              call sym_g_red_k(tt,tmp1,k1,ig)
	          call zone1_number(pnt_c(1,k1),rb0_c,ndiv_c,k11)
	          tm(k11,:,:)=tt
	        enddo
	      enddo   !! over k
	      call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nfun_red**2,tm,1)
          do ir=1,nqdiv_c
	        g_mm(:,:,ir,it,ind_tau)=tm(ir,:,:)/nqdiv_c
	      enddo
        enddo  !! over it
      enddo   !! over ind_tau
	  deallocate(g,tm,tmp,tmp1,tt)
      end
