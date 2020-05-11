      subroutine dp_bnd_uniform(k,g_omega,lambda_dyn,pol,
     &                          lambda_stat,g_x,ex,ispin)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin,k
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c),
     &                      ex(nbndf_bnd,npnt_c)
      complex*16, intent(in) :: lambda_dyn(nbndf_bnd,nbndf_bnd,
     &                                     0:n_omega),
     &                          lambda_stat(nbndf_bnd,nbndf_bnd),
     &                          g_omega(nbndf_bnd,nbndf_bnd,0:n_omega,
     &                                  npnt_c)
	  complex*16, intent(inout) :: pol
      integer :: i_omega,nn,n,k0,nk,i,j,kf0
      real*8 :: om,de
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
      allocate(tmp(nk,nk))
      allocate(tmp1(nk,nk))
      allocate(kx(nk,nk))
      allocate(vx(nk,nk)) 
      allocate(tt(nk,nk))
      tt=(0.d0,0.d0)      
c ------- Asymptotic preparations -------------------------------------
      allocate(vv(nk,nk))
      vv=lambda_stat(1:nk,1:nk)
	  do i_omega=0,n_omega
	    om=w_omega(i_omega)
        vx=vv+lambda_dyn(1:nk,1:nk,i_omega)
	    call zgemm('n','n',nk,nk,nk,(1.d0,0.d0),vx,nk,
     &             g_omega(1,1,i_omega,k0),n,(0.d0,0.d0),tmp,nk)
	    call zgemm('n','n',nk,nk,nk,(1.d0,0.d0),
     &	           g_omega(1,1,i_omega,k0),n,tmp,nk,(0.d0,0.d0),
     &             tmp1,nk)
c -------- Exchange part -------------------------
        do j=1,nk
          do i=1,nk
            kx(i,j)=vv(i,j)/dcmplx(-ex(i,k0),om)
     &                     /dcmplx(-ex(j,k0),om)
          enddo
        enddo
        tmp1=tmp1-kx
        tmp1=tmp1+conjg(transpose(tmp1))
        tt=tt+tmp1*tr(i_omega)
      enddo  !! over i_omega
c --------- Add the asymptotic correction --------------------------
      do j=1,nk
        do i=1,nk
          de=ex(j,k0)-ex(i,k0)
          if(abs(de)<1.d-8) then
            kx(i,j)=vv(i,j)*g_x(i,n_tau,k0)
     &                     *betta_t*(1.d0+g_x(i,n_tau,k0))
          else
            kx(i,j)=vv(i,j)*(g_x(i,n_tau,k0)-g_x(j,n_tau,k0))/de
          endif
        enddo
      enddo
      tt=tt+kx
      do i=1,nk
        pol=pol+tt(i,i)/nqdiv_c/sqrt(amega)
      enddo
      deallocate(vv,tmp,tmp1,kx,vx,tr,tt)
      end
