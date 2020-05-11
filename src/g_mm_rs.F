      subroutine g_mm_rs(g_tau,g_mm,z_red,ispin)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin
	  real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c)
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd,nqdiv_c)
      real*8, intent(out) :: g_mm(nrel,nfun_red,nfun_red,0:n_tau,
     &                            nqdiv_c)
      integer :: k,ind_tau,it,i_tau,nb,k0,k1,k11,ir,ii,kf,ka
      complex*16, allocatable :: g(:,:),tm(:,:,:),tmp(:,:),tmp1(:,:)
	  allocate(g(nbndf_bnd,nbndf_bnd))
	  allocate(tm(nqdiv_c,nfun_red,nfun_red))
	  allocate(tmp(nfun_red,nfun_red))
	  allocate(tmp1(nbndf_bnd,nfun_red))
	  g_mm=0.d0
	  do ind_tau=1,ndim3_tau
        i_tau=ndim3_tau*me_t+ind_tau-1
        do it=1,2
          if(it==2) i_tau=n_tau-i_tau
          tm=(0.d0,0.d0)
	      do k=1,npnt_c
            kf=k_a_from_c(k)
	        nb=n_low_bnd(kf,ispin)
c -------------------- G on K_TAU ------------------------------------
            call ferm_unpack_tau(g,g_tau(1,1,1,ind_tau,k),nb,nbndf_bnd,
     &                           nbndf_bnd,it)
            do k0=1,k_star_c(k)
	          k1=k_list_c(k0,k)
              ka=k_a_from_c(k1)
	          call zone1_number(pnt_c(1,k1),rb0_c,ndiv_c,k11)
              call zgemm('n','c',nb,nfun_red,nb,(1.d0,0.d0),g,nbndf_bnd,
     &                   z_red(1,1,k1),nfun_red,(0.d0,0.d0),tmp1,
     &                   nbndf_bnd)
	          call zgemm('n','n',nfun_red,nfun_red,nb,(1.d0,0.d0),
     &                   z_red(1,1,k1),nfun_red,tmp1,nbndf_bnd,
     &                   (0.d0,0.d0),tmp,nfun_red)
	          tm(k11,:,:)=tmp
	        enddo
	      enddo   !! over k
	      call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nfun_red**2,tm,1)
          do ir=1,nqdiv_c
	        g_mm(1,:,:,i_tau,ir)=real(tm(ir,:,:))/nqdiv_c
            if(irel==2) then
	          g_mm(2,:,:,i_tau,ir)=imag(tm(ir,:,:))/nqdiv_c
            endif
	      enddo
        enddo  !! over it
      enddo   !! over ind_tau
	  deallocate(g,tm,tmp,tmp1)
	  if(nproc_t/=1) then
        ii=nqdiv_c*(n_tau+1)
	    call dgop(g_mm,nrel*nfun_red**2*ii,'  +',comm_t)
      endif
      end
