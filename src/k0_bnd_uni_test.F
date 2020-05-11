      subroutine k0_bnd_uni_test(k0_pw,g_omega,g_x,ex)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c,nspin),
     &                      ex(nbndf_bnd,npnt_c,nspin),
     &                      g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,
     &                              npnt_c,nspin)
      real*8, intent(out) :: k0_pw(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                             npnt_c,nspin)
      integer :: i_omega,k,ispin,i_tau,nk,i,ka,ind_omega,ind_tau,it
      real*8 :: om,x1,x2,cf(2),ta
      complex*16 :: kx
      real*8, allocatable :: tmp1(:,:,:,:)
      complex*16, allocatable :: tmp(:,:),tmp2(:,:)
      k0_pw=0.d0
      do ispin=1,nspin
        do k=1,npnt_c
          ka=k_a_from_c(k)
          nk=n_low_bnd(ka,ispin)
          allocate(tmp(nk,nk))
          allocate(tmp2(nk,nk))
          allocate(tmp1(nk,nk,2,ndim3_omega))
          tmp1=0.d0
          do ind_omega=1,ndim3_omega
            i_omega=me_t*ndim3_omega+ind_omega-1
            om=w_omega(i_omega)
            call ferm_unpack_omega(tmp,g_omega(1,1,1,ind_omega,k,ispin),
     &                             nk,nbndf_bnd,nk)
            call zgemm('n','n',nk,nk,nk,(-1.d0,0.d0),tmp,nk,tmp,nk,
     &                 (0.d0,0.d0),tmp2,nk)
c -------- Exchange part -------------------------
            do i=1,nk
              kx=-(1.d0,0.d0)/dcmplx(-ex(i,k,ispin),om)**2
              tmp2(i,i)=tmp2(i,i)-kx
            enddo
            call ferm_pack_omega(tmp2,tmp1(1,1,1,ind_omega),nk,nk,nk)
          enddo  !! over ind_omega
          call omega_to_tau_baa(tmp1,nk,k0_pw(1,1,1,1,k,ispin),
     &                          nbndf_bnd,nk)
          deallocate(tmp,tmp1,tmp2)
c --------- Add the asymptotic correction --------------------------
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            do i=1,nk
              ta=tau_mesh(i_tau)
              x1=g_x(i,i_tau,k,ispin)*(ta+betta_t*g_x(i,n_tau,k,ispin))
              x2=g_x(i,n_tau-i_tau,k,ispin)
     &                   *(betta_t-ta+betta_t*g_x(i,n_tau,k,ispin))
              cf(1)=x1-x2
              cf(2)=x1+x2
              do it=1,2
                k0_pw(i,i,it,ind_tau,k,ispin)=
     &              k0_pw(i,i,it,ind_tau,k,ispin)+cf(it)
              enddo
            enddo
          enddo  !! over ind_tau
        enddo   !! over k
      enddo   !! over ispin
      k0_pw=k0_pw/sqrt(amega)
      end
