      subroutine k_bnd_stat(ispin,nom,ind_nu,k,iq,k_pw,lambda_stat,
     &                      g_om_nu,g_x,ex)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ispin,nom,ind_nu,k,iq
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c),
     &                      ex(nbndf_bnd,npnt_c)
      complex*16, intent(in) :: lambda_stat(nbndf_bnd,nbndf_bnd),
     &                          g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,
     &                                  npnt_c)
	complex*16, intent(out) :: k_pw(nbndf_bnd,nbndf_bnd)
      integer :: i_omega,kq,n,nn,in,k0,nk,kq0,nkq,j,i,i_nu,kf0,kfq0
      real*8 :: v(3),omnu,om,de
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),vx(:,:),vv(:,:),
     &                           kx(:,:)
	n=nbndf_bnd
	nn=n*n
	i_nu=me_t*ndim3_nu+ind_nu-1
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
      allocate(vx(nk,nkq))
      allocate(tmp(nk,nkq))
      allocate(tmp1(n,n))   
c ------- Asymptotic preparations -------------------------------------
      allocate(vv(nk,nkq))
      vv=lambda_stat(1:nk,1:nkq)
      allocate(kx(nk,nkq))
c -------------------------------------------------------------------
	do in=1,2
	  do i_omega=1,nom
	    om=w_om_adapt_nu(i_omega,i_nu)
          omnu=om-w_nu(i_nu)
          vx=vv
          tmp1=(0.d0,0.d0)
	    if(in==1) then
	      call zgemm('c','n',nk,nkq,nk,(1.d0,0.d0),
     &	             g_om_nu(1,1,2,i_omega,k0),n,vx,nk,(0.d0,0.d0),
     &	             tmp,nk)
	      call zgemm('n','c',nk,nkq,nkq,(1.d0,0.d0),tmp,nk,
     &	             g_om_nu(1,1,1,i_omega,kq0),n,(0.d0,0.d0),
     &	             tmp1,n)
c -------- Exchange part -------------------------
            do j=1,nkq
              do i=1,nk
                kx(i,j)=-vv(i,j)/dcmplx(-ex(i,k0),-omnu)
     &                          /dcmplx(-ex(j,kq0),-om)
              enddo
            enddo
	    else if(in==2) then
	      call zgemm('n','n',nk,nkq,nk,(1.d0,0.d0),
     &	             g_om_nu(1,1,1,i_omega,k0),n,vx,nk,(0.d0,0.d0),
     &	             tmp,nk)
	      call zgemm('n','n',nk,nkq,nkq,(1.d0,0.d0),tmp,nk,
     &	             g_om_nu(1,1,2,i_omega,kq0),n,(0.d0,0.d0),
     &	             tmp1,n)
c -------- Exchange part -------------------------
            do j=1,nkq
              do i=1,nk
                kx(i,j)=-vv(i,j)/dcmplx(-ex(i,k0),om)
     &                          /dcmplx(-ex(j,kq0),omnu)
              enddo
            enddo
          endif
          tmp1=-tmp1
c -------- Subtract Exchange part -------------------------
          tmp1(1:nk,1:nkq)=tmp1(1:nk,1:nkq)-kx
          k_pw=k_pw+tau_nu_from_omega_nu(i_omega,0,ind_nu)*tmp1
        enddo  !! over i_omega
c --------- Add the asymptotic correction --------------------------
        do j=1,nkq
          do i=1,nk
            de=ex(j,kq0)-ex(i,k0)
            if(i_nu==0.and.abs(de)<1.d-8) then
              kx(i,j)=0.5d0*g_x(i,0,k0)*betta_t*g_x(i,n_tau,k0)
            else if(in==1) then
              kx(i,j)=g_x(i,0,k0)/dcmplx(de,w_nu(i_nu))
            else if(in==2) then
              kx(i,j)=-g_x(j,0,kq0)/dcmplx(de,w_nu(i_nu))
            endif
            kx(i,j)=kx(i,j)*vv(i,j)
          enddo
        enddo
	  k_pw(1:nk,1:nkq)=k_pw(1:nk,1:nkq)+kx
      enddo
      deallocate(tmp,tmp1,vx,vv,kx)
      end