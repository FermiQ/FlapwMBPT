      subroutine vertex_g2_t2(iq,i_nu,s2,nom,k_pw,z_red,lambda_dyn,ch,
     &                        t_rr)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      integer, intent(in) :: iq,i_nu,nom,s2
      real*8, intent(in) :: t_rr(n_lem2_red,n_lem2_red,2,ndim3_tau,
     &                           ndimc_kk(me_k+1),nspin,nspin)
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,2,ndim3_tau,2,
     &                               nqdiv_c,nspin),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(out) :: lambda_dyn(nbndf_bnd,nbndf_bnd,nom,2,
     &                                      ndimc_kk(me_k+1),nspin)
      integer :: ispin,it,k,il,ind_tau,ind_r,jspin
      real*8, allocatable :: tmp(:,:,:)
      complex*16, allocatable :: lambda_tau(:,:,:,:,:,:,:),k_rs(:,:,:),
     &                           kcom(:,:,:)
      allocate(lambda_tau(nbndf_bnd,nbndf_bnd,ndimc_kk(me_k+1),2,2,
     &                    ndim3_tau,nspin))
      lambda_tau=(0.d0,0.d0)
      allocate(tmp(n_lem2_red,n_lem2_red,ndimc_kk(me_k+1)))
      allocate(kcom(nbndf_bnd,nbndf_bnd,nqdiv_c))
      allocate(k_rs(nfun_red,nfun_red,ndimc_kk(me_k+1)))
      do il=1,2
        do ind_tau=1,ndim3_tau
          do jspin=1,nspin
            do it=1,2
              do k=1,nqdiv_c
                kcom(:,:,k)=k_pw(:,:,it,ind_tau,il,k,jspin)
              enddo
              call k_real_space_sf(iq,kcom,k_rs,z_red,jspin)
              do ispin=1,nspin
                do ind_r=1,ndimc_kk(me_k+1)
                  tmp(:,:,ind_r)=t_rr(:,:,it,ind_tau,ind_r,jspin,ispin)
                enddo
                call vrt_g2_t2_mm(ispin,iq,k_rs,z_red,tmp,
     &                            lambda_tau(1,1,1,il,it,ind_tau,ispin),
     &                            ch)
              enddo
            enddo
          enddo
        enddo
      enddo  !! over ind_tau
      deallocate(tmp,kcom,k_rs)
      do ispin=1,nspin
        call vertex_tau_to_omega_nu_par(lambda_tau(1,1,1,1,1,1,ispin),
     &                                  lambda_dyn(1,1,1,1,1,ispin),
     &                                  i_nu,iq,nom,ispin)
      enddo
      deallocate(lambda_tau)
      end
