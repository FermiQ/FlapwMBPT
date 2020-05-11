      subroutine k_bnd(ispin,nom,ind_nu,k,iq,lambda_dyn,k_pw,
     &                 lambda_stat,g_om_nu,g_x,ex)
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
      complex*16, intent(in) :: lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2),
     &                          lambda_stat(nbndf_bnd,nbndf_bnd),
     &                          g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,
     &                                  npnt_c)
	  complex*16, intent(out) :: k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,2)
      integer :: i_omega,kq,n,nn,in,i_tau,k0,nk,kq0,nkq,j,i,i_nu,kf0,
     &           kfq0
      real*8 :: v(3),omnu,om,de
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:),vx(:,:),vv(:,:),
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
      allocate(tmp1(n,n,nom))   
c ------- Asymptotic preparations -------------------------------------
      allocate(vv(nk,nkq))
      vv=lambda_stat(1:nk,1:nkq)
      allocate(kx(nk,nkq))
c -------------------------------------------------------------------
	  do in=1,2
        tmp1=(0.d0,0.d0)
	    do i_omega=1,nom
	      om=w_om_adapt_nu(i_omega,i_nu)
          omnu=om-w_nu(i_nu)
          vx=vv+lambda_dyn(1:nk,1:nkq,i_omega,in)
	      if(in==1) then
	        call zgemm('c','n',nk,nkq,nk,(1.d0,0.d0),
     &	               g_om_nu(1,1,2,i_omega,k0),n,vx,nk,(0.d0,0.d0),
     &	               tmp,nk)
	        call zgemm('n','c',nk,nkq,nkq,(1.d0,0.d0),tmp,nk,
     &	               g_om_nu(1,1,1,i_omega,kq0),n,(0.d0,0.d0),
     &	               tmp1(1,1,i_omega),n)
c -------- Exchange part -------------------------
            do j=1,nkq
              do i=1,nk
                kx(i,j)=-vv(i,j)/dcmplx(-ex(i,k0),-omnu)
     &                          /dcmplx(-ex(j,kq0),-om)
              enddo
            enddo
	      else if(in==2) then
	        call zgemm('n','n',nk,nkq,nk,(1.d0,0.d0),
     &	               g_om_nu(1,1,1,i_omega,k0),n,vx,nk,(0.d0,0.d0),
     &	               tmp,nk)
	        call zgemm('n','n',nk,nkq,nkq,(1.d0,0.d0),tmp,nk,
     &	               g_om_nu(1,1,2,i_omega,kq0),n,(0.d0,0.d0),
     &	               tmp1(1,1,i_omega),n)
c -------- Exchange part -------------------------
            do j=1,nkq
              do i=1,nk
                kx(i,j)=-vv(i,j)/dcmplx(-ex(i,k0),om)
     &                          /dcmplx(-ex(j,kq0),omnu)
              enddo
            enddo
          endif
          tmp1(:,:,i_omega)=-tmp1(:,:,i_omega)
c -------- Subtract Exchange part -------------------------
          tmp1(1:nk,1:nkq,i_omega)=tmp1(1:nk,1:nkq,i_omega)-kx
          if(in==1) tmp1(:,:,i_omega)=conjg(tmp1(:,:,i_omega))
        enddo  !! over i_omega
        call from_omega_nu_to_tau_nu_ab(nn,ind_nu,tmp1,k_pw(1,1,0,in),
     &                                  in)
c --------- Add the asymptotic correction --------------------------
        do i_tau=0,n_tau
          do j=1,nkq
            do i=1,nk
              de=ex(j,kq0)-ex(i,k0)
              if(i_nu==0.and.abs(de)<1.d-8) then
                kx(i,j)=0.5d0*g_x(i,i_tau,k0)
     &            *(betta_t*g_x(i,n_tau,k0)+tau_mesh(i_tau))
              else if(in==1) then
                kx(i,j)=g_x(i,i_tau,k0)
     &                 /dcmplx(de,w_nu(i_nu))
              else if(in==2) then
                kx(i,j)=-g_x(j,i_tau,kq0)
     &                  /dcmplx(de,w_nu(i_nu))
              endif
              kx(i,j)=kx(i,j)*vv(i,j)
            enddo
          enddo
	      k_pw(1:nk,1:nkq,i_tau,in)=k_pw(1:nk,1:nkq,i_tau,in)+kx
        enddo
      enddo
      deallocate(tmp,tmp1,vx,vv,kx)
      end
