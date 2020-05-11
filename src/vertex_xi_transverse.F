      subroutine vertex_xi_transverse(n,nl2,nl2m,ind_nl2,v4_q,w4_q_tau,
     &                                k_pw,iq,nom,lambda_dyn,
     &                                lambda_stat,ind_nu,isp1,isp2,key)
c     key =0  for P_00
c     key =1  for P_zz
c     key =2  for P_+-
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: nl2,n,iq,ind_nl2(n,n,nqdiv),nom,ind_nu,
     &                       isp1,isp2,key,nl2m
      complex*16, intent(in) :: v4_q(nl2m,nl2m),
     &                          w4_q_tau(nl2m,nl2m,0:n_tau/2),
     &                          k_pw(n,n,0:n_tau,2,nqdiv)
      complex*16, intent(inout) :: lambda_dyn(n,n,n_omega_b_max,2,
     &                                        nqdiv),
     &                             lambda_stat(n,n,nqdiv)
      integer :: k,kq,j,i,ii,i_tau,i1_tau,i0_tau,in
      real*8 :: v(3),fc
      complex*16, allocatable :: vx_dyn(:,:,:),vx(:),tmp(:,:),dk(:),
     &                           kt(:,:,:),lam(:,:,:)
      fc=1.d0/nqdiv
      if(key==2) fc=2.d0/nqdiv
	allocate(vx_dyn(nl2,0:n_tau,2))
	allocate(vx(nl2))
	allocate(tmp(nl2,2))
	allocate(kt(nl2,2,0:n_tau))
	allocate(dk(nl2))
      do k=1,nqdiv
	  v=pnt(:,k)-pnt(:,iq)
	  call zone1_number(v,rb0,ndiv,kq)
	  kq=index_k1(kq)     !! for K-Q
        do j=1,n_low_bnd(i_kref(kq),isp2)
          do i=1,n_low_bnd(i_kref(k),isp1)
            ii=ind_nl2(i,j,k)
            dk(ii)=k_pw(i,j,0,1,k)+k_pw(i,j,0,2,k)
            kt(ii,1,:)=k_pw(i,j,:,1,k)
            kt(ii,2,:)=k_pw(i,j,:,2,k)
          enddo
        enddo
      enddo   !! over k
	call zgemm('n','n',nl2,1,nl2,(1.d0,0.d0),v4_q,nl2,dk,nl2,
     &	         (0.d0,0.d0),vx,nl2) 
      vx=fc*vx
      do i_tau=0,n_tau
        i1_tau=n_tau-i_tau
        i0_tau=min(i_tau,i1_tau)
c ----- Dynamical part : Lambda = W*K ----------------------------
	  call zgemm('n','n',nl2,2,nl2,(1.d0,0.d0),w4_q_tau(1,1,i0_tau),
     &	           nl2,kt(1,1,i_tau),nl2,(0.d0,0.d0),tmp,nl2)
        vx_dyn(:,i_tau,1)=fc*tmp(:,1)
        vx_dyn(:,i_tau,2)=fc*tmp(:,2)
      enddo  
c -------- Transform Vrt(tau;nu) ---- > Vrt(omega;nu) ----------------
	allocate(lam(nl2,0:nom,2))
      do in=1,2
        call from_tau_nu_to_omega_nu_d(nl2,ind_nu,nom,lam(1,1,in),
     &                                 vx_dyn,in)
      enddo
      do k=1,nqdiv
	  v=pnt(:,k)-pnt(:,iq)
	  call zone1_number(v,rb0,ndiv,kq)
	  kq=index_k1(kq)     !! for K-Q
        do j=1,n_low_bnd(i_kref(kq),isp2)
          do i=1,n_low_bnd(i_kref(k),isp1)
            ii=ind_nl2(i,j,k)
            lambda_stat(i,j,k)=vx(ii)
            do in=1,2
              lambda_dyn(i,j,1:nom,in,k)=lam(ii,1:nom,in)
            enddo
          enddo
        enddo
      enddo   !! over k
      deallocate(vx_dyn,vx,dk,kt,lam,tmp)
      end