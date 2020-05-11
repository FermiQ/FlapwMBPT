      subroutine k_bnd_uni_test(k_pw,g_omega,lambda_dyn,lambda_stat,g_x,
     &                          ex)
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
     &                              npnt_c,nspin),
     &                      lambda_dyn(nbndf_bnd,nbndf_bnd,2,
     &                                 ndim3_omega,npnt_c,nspin)
      complex*16, intent(in) :: lambda_stat(nbndf_bnd,nbndf_bnd,npnt_c,
     &                                      nspin)
      real*8, intent(out) :: k_pw(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c,nspin)
      integer :: i_omega,k,ispin,i_tau,nk,i,ka,ind_omega,ind_tau,it,j
      real*8 :: om,de
      complex*16 :: kx
      real*8, allocatable :: tmp1(:,:,:,:)
      complex*16, allocatable :: tmp(:,:),tmp2(:,:),tmv(:,:),kx2(:,:,:)
      k_pw=0.d0
      do ispin=1,nspin
        do k=1,npnt_c
          ka=k_a_from_c(k)
          nk=n_low_bnd(ka,ispin)
          allocate(tmp(nk,nk))
          allocate(tmv(nk,nk))
          allocate(tmp2(nk,nk))
          allocate(tmp1(nk,nk,2,ndim3_omega))
          tmp1=(0.d0,0.d0)
          do ind_omega=1,ndim3_omega
            i_omega=me_t*ndim3_omega+ind_omega-1
            om=w_omega(i_omega)
            call ferm_unpack_omega(tmp,g_omega(1,1,1,ind_omega,k,ispin),
     &                             nk,nbndf_bnd,nk)
            call ferm_unpack_omega(tmv,
     &                             lambda_dyn(1,1,1,ind_omega,k,ispin),
     &                             nk,nbndf_bnd,nk)
            tmv=tmv+lambda_stat(1:nk,1:nk,k,ispin)
            call zgemm('n','n',nk,nk,nk,(-1.d0,0.d0),tmp,nk,tmv,nk,
     &                 (0.d0,0.d0),tmp2,nk)
            call zgemm('n','n',nk,nk,nk,(1.d0,0.d0),tmp2,nk,tmp,nk,
     &                 (0.d0,0.d0),tmv,nk)
c -------- Exchange part -------------------------
            do j=1,nk
              do i=1,nk
                kx=-lambda_stat(i,j,k,ispin)/dcmplx(-ex(i,k,ispin),om)
     &                                      /dcmplx(-ex(j,k,ispin),om)
                tmv(i,j)=tmv(i,j)-kx
              enddo
            enddo
            call ferm_pack_omega(tmv,tmp1(1,1,1,ind_omega),nk,nk,nk)
          enddo  !! over i_omega
          call omega_to_tau_baa(tmp1,nk,k_pw(1,1,1,1,k,ispin),
     &                          nbndf_bnd,nk)
          deallocate(tmp,tmp1,tmp2,tmv)
c --------- Add the asymptotic correction --------------------------
          allocate(kx2(nk,nk,2))
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              do j=1,nk
                do i=1,nk
                  de=ex(j,k,ispin)-ex(i,k,ispin)
                  if(abs(de)<1.d-8) then
                    kx=g_x(i,i_tau,k,ispin)
     &                *(betta_t*g_x(i,n_tau,k,ispin)+tau_mesh(i_tau))
                  else
                    kx=(g_x(i,i_tau,k,ispin)-g_x(j,i_tau,k,ispin))/de
                  endif
                  kx2(i,j,it)=kx*lambda_stat(i,j,k,ispin)
                enddo
              enddo
            enddo
            call ferm_pack_tau(kx2,k_pw(1,1,1,ind_tau,k,ispin),nk,nk,
     &                         nbndf_bnd)
          enddo  !! over ind_tau
          deallocate(kx2)
        enddo   !! over k
      enddo   !! over ispin
      end
