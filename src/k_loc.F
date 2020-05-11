      subroutine k_loc(ispin,nom,ind_nu,lambda_dyn,k_pw,
     &                 lambda_stat,g_om_nu,n)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ispin,nom,ind_nu,n
      complex*16, intent(in) :: lambda_dyn(n,n,nom,2),
     &                          lambda_stat(n,n),
     &                          g_om_nu(n,n,2,nom)
      complex*16, intent(out) :: k_pw(n,n,0:n_tau,2)
      integer :: i_omega,nn,in,i_nu
      real*8 :: omnu,om
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:),vx(:,:),vv(:,:)
      nn=n*n
      i_nu=me_t*ndim3_nu+ind_nu-1
c ------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
      allocate(vx(n,n))
      allocate(tmp(n,n))
      allocate(tmp1(n,n,nom))
c ------- Asymptotic preparations -------------------------------------
      allocate(vv(n,n))
      vv=lambda_stat
c -------------------------------------------------------------------
      do in=1,2
        tmp1=(0.d0,0.d0)
        do i_omega=1,nom
          om=w_om_adapt_nu(i_omega,i_nu)
          omnu=om-w_nu(i_nu)
          vx=vv+lambda_dyn(:,:,i_omega,in)
          if(in==1) then
            call zgemm('c','n',n,n,n,(1.d0,0.d0),g_om_nu(1,1,2,i_omega),
     &                 n,vx,n,(0.d0,0.d0),tmp,n)
            call zgemm('n','c',n,n,n,(1.d0,0.d0),tmp,n,
     &                 g_om_nu(1,1,1,i_omega),n,(0.d0,0.d0),
     &                 tmp1(1,1,i_omega),n)
          else if(in==2) then
            call zgemm('n','n',n,n,n,(1.d0,0.d0),
     &                 g_om_nu(1,1,1,i_omega),n,vx,n,(0.d0,0.d0),
     &                 tmp,n)
            call zgemm('n','n',n,n,n,(1.d0,0.d0),tmp,n,
     &                 g_om_nu(1,1,2,i_omega),n,(0.d0,0.d0),
     &                 tmp1(1,1,i_omega),n)
          endif
          tmp1(:,:,i_omega)=-tmp1(:,:,i_omega)
c -------- Subtract Exchange part -------------------------
          if(in==1) tmp1(:,:,i_omega)=conjg(tmp1(:,:,i_omega))
        enddo  !! over i_omega
        call from_omega_nu_to_tau_nu_ab(nn,ind_nu,tmp1,k_pw(1,1,0,in),
     &                                  in)
      enddo
      deallocate(tmp,tmp1,vx,vv)
      end
