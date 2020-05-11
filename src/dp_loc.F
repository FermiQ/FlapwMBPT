      subroutine dp_loc(nom,ind_nu,lambda_dyn,lambda_stat,
     &                  g_om_nu,fif,pol,ispin,trn,key,s2,n,np)
c     key = 0 for full P = P0 + dP
c     key = 1 for dP only
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nom,ind_nu,ispin,key,s2,n,np
      real*8, intent(in) :: fif(nrel,n,n,np)
      complex*16, intent(in) :: lambda_dyn(n,n,nom,2),
     &                          lambda_stat(n,n),
     &                          g_om_nu(n,n,2,nom),
     &                          trn(n_omega1_max)
      complex*16, intent(inout) :: pol(np)
      integer :: i_omega,nn,in,ii,i,j,i_nu
      real*8 :: om,omnu
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),vx(:,:),tt(:,:),
     &                           xx(:,:),ff(:,:,:)
      i_nu=me_t*ndim3_nu+ind_nu-1
      nn=n*n
c ------------------------------------------------------------------
      allocate(tmp(n,n))
      allocate(tmp1(n,n))
      allocate(vx(n,n))
      allocate(ff(n,n,np))
      if(irel/=2) ff=fif(1,:,:,:)
      if(irel==2) ff=dcmplx(fif(1,:,:,:),fif(2,:,:,:))
c ------- Asymptotic preparations -------------------------------------
      allocate(xx(n,n))
      xx=lambda_stat
      if(key==0) xx=xx+ff(:,:,s2)
c -------------------------------------------------------------------
      allocate(tt(n,n))
      tt=(0.d0,0.d0)
      do in=1,2
        do i_omega=1,nom
          om=w_om_adapt_nu(i_omega,i_nu)
          omnu=om-w_nu(i_nu)
          vx=xx+lambda_dyn(:,:,i_omega,in)
          if(in==1) then
            call zgemm('c','n',n,n,n,(1.d0,0.d0),
     &                 g_om_nu(1,1,2,i_omega),n,vx,n,(0.d0,0.d0),
     &                 tmp,n)
            call zgemm('n','c',n,n,n,(1.d0,0.d0),tmp,n,
     &                 g_om_nu(1,1,1,i_omega),n,(0.d0,0.d0),tmp1,
     &                 n)
c -------------------------------------------------------------------
            tt=tt+tmp1*conjg(trn(i_omega))
          else if(in==2) then
            call zgemm('n','n',n,n,n,(1.d0,0.d0),
     &                 g_om_nu(1,1,1,i_omega),n,vx,n,(0.d0,0.d0),
     &                 tmp,n)
            call zgemm('n','n',n,n,n,(1.d0,0.d0),tmp,n,
     &                 g_om_nu(1,1,2,i_omega),n,(0.d0,0.d0),tmp1,n)
c -------------------------------------------------------------------
            tt=tt+tmp1*trn(i_omega)
          endif
        enddo  !! over i_omega
      enddo
      do ii=1,np
        do j=1,n
          do i=1,n
            pol(ii)=pol(ii)+conjg(ff(i,j,ii))*tt(i,j)
          enddo
        enddo
      enddo
      deallocate(tmp,tmp1,vx,tt,xx)
      end
