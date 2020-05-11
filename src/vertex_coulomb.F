      subroutine vertex_coulomb(iq,np,k_pw,ps3,v2_q,lambda_stat)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,np
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2,
     &                               nqdiv_c,nspin),
     &                          ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,
     &                              nqdiv_c,nspin),
     &                          v2_q(n_pbtot_red,n_pbtot_red)
      complex*16, intent(inout) :: lambda_stat(nbndf_bnd,nbndf_bnd,
     &                                         nqdiv_c,nspin)
      integer :: k0,nk,kq,kq0,i,j,k,nkq,ispin,kf0,kfq0,ii
      real*8 :: v(3)
      complex*16, allocatable :: vp(:),tt(:,:),ps(:)
      allocate(tt(nbndf_bnd,nbndf_bnd))
      allocate(ps(np))
      ps=(0.d0,0.d0)
	  do ispin=1,nspin
        do k=1,nqdiv_c
	      k0=i_kref_c(k)
          kf0=k_a_from_c(k0)
	      nk=n_low_bnd(kf0,ispin)
	      v=pnt_c(:,k)-pnt_c(:,iq)
	      call zone1_number(v,rb0_c,ndiv_c,kq)
	      kq=index_k1_c(kq)     !! for K-Q
	      kq0=i_kref_c(kq)
          kfq0=k_a_from_c(kq0)
	      nkq=n_low_bnd(kfq0,ispin)
	      tt(1:nk,1:nkq)=k_pw(1:nk,1:nkq,0,1,k,ispin)
     &                  +k_pw(1:nk,1:nkq,0,2,k,ispin)
          do ii=1,np
            do j=1,nkq
              do i=1,nk
                ps(ii)=ps(ii)-conjg(ps3(i,j,ii,k,ispin))*tt(i,j)
              enddo
            enddo
          enddo
        enddo  !! over k
      enddo  !! over ispin
      if(nspin==1.and.nrel==1) ps=2.d0*ps
      ps=ps/nqdiv_c
      allocate(vp(np))
      vp=(0.d0,0.d0)
      do j=1,np
        do i=1,np
          vp(i)=vp(i)+v2_q(i,j)*ps(j)
        enddo
      enddo
	  do ispin=1,nspin
        do k=1,nqdiv_c
	      k0=i_kref_c(k)
          kf0=k_a_from_c(k0)
	      nk=n_low_bnd(kf0,ispin)
	      v=pnt_c(:,k)-pnt_c(:,iq)
	      call zone1_number(v,rb0_c,ndiv_c,kq)
	      kq=index_k1_c(kq)     !! for K-Q
	      kq0=i_kref_c(kq)
          kfq0=k_a_from_c(kq0)
	      nkq=n_low_bnd(kfq0,ispin)
          do ii=1,np
            do j=1,nkq
              do i=1,nk
                lambda_stat(i,j,k,ispin)=lambda_stat(i,j,k,ispin)
     &                                  +ps3(i,j,ii,k,ispin)*vp(ii)
              enddo
            enddo
          enddo
        enddo  !! over k
      enddo  !! over ispin
      deallocate(vp,tt,ps)
      end
      
      
      
      
      subroutine p0_red(iq,np,ps3,g_tau,xi0)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: iq,np
	real*8, intent(out) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                             npnt_c,nspin)
      complex*16, intent(in) :: ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,
     &                              nqdiv_c,nspin)
      real*8, intent(out) :: xi0(np,np,ndim3_nu)
      integer :: k0,nk,kq,kq0,k,nkq,ispin,kf0,kfq0,ii,ind_tau
      real*8 :: v(3)
      real*8, allocatable :: xit(:,:,:)
      complex*16, allocatable :: gt(:,:),gtq(:,:),tt(:,:),tt1(:,:,:),
     &                           tt0(:,:,:),tmp(:,:)
      allocate(xit(np,np,ndim3_tau))
      allocate(gt(nbndf_bnd,nbndf_bnd))
      allocate(gtq(nbndf_bnd,nbndf_bnd))
      allocate(tt(nbndf_bnd,nbndf_bnd))
      allocate(tmp(np,np))
      xit=0.d0
	do ispin=1,nspin
        do k=1,nqdiv_c
	    k0=i_kref_c(k)
          kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
	    v=pnt_c(:,k)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K-Q
	    kq0=i_kref_c(kq)
          kfq0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kfq0,ispin)
          allocate(tt0(nk,nkq,np))
          allocate(tt1(nk,nkq,np))
          tt0=ps3(1:nk,1:nkq,1:np,k,ispin)
          do ind_tau=1,ndim3_tau
            call ferm_unpack_tau(gt,g_tau(1,1,1,ind_tau,k0,ispin),nk,
     &                           nbndf_bnd,nbndf_bnd,1)
            call ferm_unpack_tau(gtq,g_tau(1,1,1,ind_tau,kq0,ispin),nkq,
     &                           nbndf_bnd,nbndf_bnd,2)
            do ii=1,np
	        call zgemm('n','n',nk,nkq,nk,(1.d0,0.d0),gt,nbndf_bnd,
     &	               tt0(1,1,ii),nk,(0.d0,0.d0),tt,nbndf_bnd)
	        call zgemm('n','n',nk,nkq,nkq,(1.d0,0.d0),tt,nbndf_bnd,
     &	               gtq,nbndf_bnd,(0.d0,0.d0),tt1(1,1,ii),nk)
            enddo
	      call zgemm('c','n',np,np,nk*nkq,(-1.d0,0.d0),tt0,nk*nkq,
     &	             tt1,nk*nkq,(0.d0,0.d0),tmp,np)
            call boson_pack_tau(1.d0,tmp,np,1.d0,xit(1,1,ind_tau),np,np)
          enddo
          deallocate(tt0,tt1)
        enddo  !! over k
      enddo  !! over ispin
      if(nspin==1.and.nrel==1) xit=2.d0*xit
      xit=xit/nqdiv_c
      call tau_to_nu_ba(xi0,np,xit,np,np)
      deallocate(gt,gtq,xit,tt,tmp)
      end
