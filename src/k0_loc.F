      subroutine k0_loc(ind_nu,k0_pw,g_om_nu,fif,s2,nom,n,np)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,s2,nom,n,np
      real*8, intent(in) :: fif(nrel,n,n,np,nspin)
      complex*16, intent(in) :: g_om_nu(n,n,2,nom,nspin)
      complex*16, intent(out) :: k0_pw(n,n,0:n_tau,2,nspin)
      integer :: i_omega,ispin,nn,in,i_nu
      real*8 :: om,omnu
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:),ff(:,:)
      i_nu=me_t*ndim3_nu+ind_nu-1
      nn=n*n
      allocate(tmp(n,n))
      allocate(tmp1(n,n,n_omega1_max))
      allocate(ff(n,n))
c ------------------------------------------------------------------
      k0_pw=(0.d0,0.d0)
      do ispin=1,nspin
        if(irel/=2) ff=fif(1,:,:,s2,ispin)
        if(irel==2) ff=dcmplx(fif(1,:,:,s2,ispin),fif(2,:,:,s2,ispin))
        do in=1,2
          tmp1=(0.d0,0.d0)
          do i_omega=1,nom
            om=w_om_adapt_nu(i_omega,i_nu)
            omnu=om-w_nu(i_nu)
            if(in==1) then
              call zgemm('n','c',n,n,n,(1.d0,0.d0),ff,n,
     &                   g_om_nu(1,1,1,i_omega,ispin),n,(0.d0,0.d0),tmp,
     &                   n)
              call zgemm('c','n',n,n,n,(1.d0,0.d0),
     &                   g_om_nu(1,1,2,i_omega,ispin),n,tmp,n,
     &                   (0.d0,0.d0),tmp1(1,1,i_omega),n)
            else if(in==2) then
              call zgemm('n','n',n,n,n,(1.d0,0.d0),ff,n,
     &                   g_om_nu(1,1,2,i_omega,ispin),n,(0.d0,0.d0),tmp,
     &                   n)
              call zgemm('n','n',n,n,n,(1.d0,0.d0),
     &                   g_om_nu(1,1,1,i_omega,ispin),n,tmp,n,
     &                   (0.d0,0.d0),tmp1(1,1,i_omega),n)
            endif
            tmp1(:,:,i_omega)=-tmp1(:,:,i_omega)
            if(in==1) tmp1(:,:,i_omega)=conjg(tmp1(:,:,i_omega))
          enddo  !! over i_omega
          call from_omega_nu_to_tau_nu_ab(nn,ind_nu,tmp1,
     &                                    k0_pw(1,1,0,in,ispin),in)
        enddo   !! over in
      enddo   !! over ispin
      deallocate(tmp,tmp1,ff)
      end
