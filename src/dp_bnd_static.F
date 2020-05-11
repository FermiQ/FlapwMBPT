      subroutine dp_bnd_static(k,iq,g_omega,lambda_dyn,ps3,pol,
     &                         lambda_stat,g_x,ex,ispin,n0)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin,k,iq,n0
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c),
     &                      ex(nbndf_bnd,npnt_c)
      complex*16, intent(in) :: lambda_dyn(nbndf_bnd,nbndf_bnd,
     &                                     0:n_omega,2),
     &                          lambda_stat(nbndf_bnd,nbndf_bnd),
     &                          g_omega(nbndf_bnd,nbndf_bnd,0:n_omega,
     &                                  npnt_c),
     &                          ps3(nbndf_bnd,nbndf_bnd,n0)
	  complex*16, intent(inout) :: pol(n0)
      integer :: i_omega,kq,nn,n,k0,nk,kq0,ii,nkq,i,j,kf0,kfq0,iqa
      real*8 :: v(3),om,de
      real*8, allocatable :: tr(:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),kx(:,:),vv(:,:),
     &                           vx(:,:),tt(:,:)
      n=nbndf_bnd
      nn=n*n
      allocate(tr(0:n_omega))
      tr=tau_from_omega(:,0,1)/4
c ------------------------------------------------------------------
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
      allocate(tmp1(nk,nkq))
      allocate(kx(nk,nkq))
      allocate(vx(nk,nkq)) 
      allocate(tt(nk,nkq))
      tt=(0.d0,0.d0)      
c ------- Asymptotic preparations -------------------------------------
      allocate(vv(nk,nkq))
      vv=lambda_stat(1:nk,1:nkq)
	  do i_omega=0,n_omega
	    om=w_omega(i_omega)
c -------------- Negative omega (in=1) --------------------------------
        vx=vv+lambda_dyn(1:nk,1:nkq,i_omega,1)
	    call zgemm('n','c',nk,nkq,nkq,(1.d0,0.d0),vx,nk,
     &             g_omega(1,1,i_omega,kq0),n,(0.d0,0.d0),tmp,nk)
	    call zgemm('c','n',nk,nkq,nk,(1.d0,0.d0),
     &	           g_omega(1,1,i_omega,k0),n,tmp,nk,(0.d0,0.d0),
     &             tmp1,nk)
c -------- Exchange part -------------------------
        do j=1,nkq
          do i=1,nk
            kx(i,j)=vv(i,j)/dcmplx(-ex(i,k0),-om)
     &                     /dcmplx(-ex(j,kq0),-om)
          enddo
        enddo
        tmp1=tmp1-kx
        tt=tt+tmp1*tr(i_omega)
c -------------- Positive omega (in=2) --------------------------------
        vx=vv+lambda_dyn(1:nk,1:nkq,i_omega,2)
	    call zgemm('n','n',nk,nkq,nkq,(1.d0,0.d0),vx,nk,
     &             g_omega(1,1,i_omega,kq0),n,(0.d0,0.d0),tmp,nk)
	    call zgemm('n','n',nk,nkq,nk,(1.d0,0.d0),
     &	           g_omega(1,1,i_omega,k0),n,tmp,nk,(0.d0,0.d0),
     &             tmp1,nk)
c -------- Exchange part -------------------------
        do j=1,nkq
          do i=1,nk
            kx(i,j)=vv(i,j)/dcmplx(-ex(i,k0),om)
     &                     /dcmplx(-ex(j,kq0),om)
          enddo
        enddo
        tmp1=tmp1-kx
        tt=tt+tmp1*tr(i_omega)
      enddo  !! over i_omega
c --------- Add the asymptotic correction --------------------------
      do j=1,nkq
        do i=1,nk
          de=ex(j,kq0)-ex(i,k0)
          if(abs(de)<1.d-8) then
            kx(i,j)=vv(i,j)*g_x(i,n_tau,k0)
     &                     *betta_t*(1.d0+g_x(i,n_tau,k0))
          else
            kx(i,j)=vv(i,j)*(g_x(i,n_tau,k0)-g_x(j,n_tau,kq0))/de
          endif
        enddo
      enddo
      tt=tt+kx
      iqa=k_a_from_c(iq)
      do ii=1,min(n0,n_pbmt_red+nplwgw_red(iqa))
        do j=1,nkq
          do i=1,nk
            pol(ii)=pol(ii)+conjg(ps3(i,j,ii))*tt(i,j)/nqdiv_c
          enddo
        enddo
      enddo
      deallocate(vv,tmp,tmp1,kx,vx,tr,tt)
      end
