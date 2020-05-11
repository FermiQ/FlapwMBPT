      subroutine k_bnd_static(ispin,k,iq,k_pw,g_omega,lambda_dyn,
     &                        lambda_stat,g_x,ex)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ispin,k,iq
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c),
     &                      ex(nbndf_bnd,npnt_c)
      complex*16, intent(in) :: lambda_dyn(nbndf_bnd,nbndf_bnd,
     &                                     0:n_omega,2),
     &                          lambda_stat(nbndf_bnd,nbndf_bnd),
     &                          g_omega(nbndf_bnd,nbndf_bnd,0:n_omega,
     &                                  npnt_c)
	complex*16, intent(out) :: k_pw(nbndf_bnd,nbndf_bnd,0:n_tau)
      integer :: i_omega,kq,nn,n,i_tau,k0,nk,kq0,nkq,i,j,kf0,kfq0
      real*8 :: v(3),om,de
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:,:),kx(:,:),
     &                           a_tau(:,:,:),b_tau(:,:,:),vv(:,:),
     &                           vx(:,:)
      n=nbndf_bnd
      nn=n*n
c ------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
	k0=i_kref_c(k)
      kf0=k_a_from_c(k0)
	nk=n_low_bnd(kf0,ispin)
	v=pnt_c(:,k)-pnt_c(:,iq)
	call zone1_number(v,rb0_c,ndiv_c,kq)
	kq=index_k1_c(kq)     !! for K-Q
	kq0=i_kref_c(kq)
      kfq0=k_a_from_c(kq0)
	nkq=n_low_bnd(kfq0,ispin)
      allocate(tmp(nk,nkq))
      allocate(tmp1(nk,nkq,0:n_omega,2))
      allocate(kx(nk,nkq))
      allocate(vx(nk,nkq))  
c ------- Asymptotic preparations -------------------------------------
      allocate(vv(nk,nkq))
      vv=lambda_stat(1:nk,1:nkq)
	tmp1=(0.d0,0.d0)    
	do i_omega=0,n_omega
	  om=w_omega(i_omega)
c -------------- Negative omega (in=1) --------------------------------
        vx=vv+lambda_dyn(1:nk,1:nkq,i_omega,1)
	  call zgemm('n','c',nk,nkq,nkq,(1.d0,0.d0),vx,nk,
     &             g_omega(1,1,i_omega,kq0),n,(0.d0,0.d0),tmp,nk)
	  call zgemm('c','n',nk,nkq,nk,(-1.d0,0.d0),
     &	         g_omega(1,1,i_omega,k0),n,tmp,nk,(0.d0,0.d0),
     &             tmp1(1,1,i_omega,1),nk)
c -------- Exchange part -------------------------
        do j=1,nkq
          do i=1,nk
            kx(i,j)=-vv(i,j)/dcmplx(-ex(i,k0),-om)
     &                      /dcmplx(-ex(j,kq0),-om)
          enddo
        enddo
        tmp1(:,:,i_omega,1)=tmp1(:,:,i_omega,1)-kx
c -------------- Positive omega (in=2) --------------------------------
        vx=vv+lambda_dyn(1:nk,1:nkq,i_omega,2)
	  call zgemm('n','n',nk,nkq,nkq,(1.d0,0.d0),vx,nk,
     &             g_omega(1,1,i_omega,kq0),n,(0.d0,0.d0),tmp,nk)
	  call zgemm('n','n',nk,nkq,nk,(-1.d0,0.d0),
     &	         g_omega(1,1,i_omega,k0),n,tmp,nk,(0.d0,0.d0),
     &             tmp1(1,1,i_omega,2),nk)
c -------- Exchange part -------------------------
        do j=1,nkq
          do i=1,nk
            kx(i,j)=-vv(i,j)/dcmplx(-ex(i,k0),om)
     &                      /dcmplx(-ex(j,kq0),om)
          enddo
        enddo
        tmp1(:,:,i_omega,2)=tmp1(:,:,i_omega,2)-kx
c --------- Transform to A and B form ---------------------------------
c --------- A=tmp1(+) + tmp1(-) ---------------------------------------
c --------- B=i{-tmp1(+) + tmp1(-)} -----------------------------------
        tmp=tmp1(:,:,i_omega,2)+tmp1(:,:,i_omega,1)
        tmp1(:,:,i_omega,1)=(0.d0,1.d0)*(-tmp1(:,:,i_omega,2)
     &                                   +tmp1(:,:,i_omega,1))
        tmp1(:,:,i_omega,2)=tmp
      enddo  !! over i_omega
c --------- Get A(tau) and B(tau) -------------------------------------
      allocate(a_tau(nk,nkq,0:n_tau/2))
      allocate(b_tau(nk,nkq,0:n_tau/2))
	call dgemm('n','n',2*nk*nkq,n_tau/2+1,n_omega+1,1.d0,
     &           tmp1(1,1,0,2),2*nk*nkq,tau_from_omega(0,0,1),
     &           n_omega+1,0.d0,a_tau,2*nk*nkq)
	call dgemm('n','n',2*nk*nkq,n_tau/2+1,n_omega+1,1.d0,
     &           tmp1(1,1,0,1),2*nk*nkq,tau_from_omega(0,0,2),
     &           n_omega+1,0.d0,b_tau,2*nk*nkq)
c -------- Get tau and b-tau components of K0_PW ----------------------
      do i_tau=0,n_tau/2
        k_pw(1:nk,1:nkq,i_tau)=
     &    0.25d0*(a_tau(:,:,i_tau)+b_tau(:,:,i_tau))
        k_pw(1:nk,1:nkq,n_tau-i_tau)=
     &    0.25d0*(b_tau(:,:,i_tau)-a_tau(:,:,i_tau))
      enddo
c --------- Add the asymptotic correction --------------------------
      do i_tau=0,n_tau
        do j=1,nkq
          do i=1,nk
            de=ex(j,kq0)-ex(i,k0)
            if(abs(de)<1.d-8) then
              kx(i,j)=g_x(i,i_tau,k0)
     &               *(betta_t*g_x(i,n_tau,k0)+tau_mesh(i_tau))
            else
              kx(i,j)=(g_x(i,i_tau,k0)-g_x(j,i_tau,kq0))/de
            endif
            kx(i,j)=kx(i,j)*vv(i,j)
          enddo
        enddo
	  k_pw(1:nk,1:nkq,i_tau)=k_pw(1:nk,1:nkq,i_tau)+kx
      enddo
      deallocate(vv,tmp,tmp1,kx,vx)
      end