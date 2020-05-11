      subroutine vertex_vk_loc(k_pw,vloc,lambda_stat,nsp,fif,n,np)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,np,nsp
      real*8, intent(in) :: vloc(n_pbmtm_red,n_pbmtm_red),
     &                      fif(nrel,n,n,np,nspin)
      complex*16, intent(in) :: k_pw(n,n,nsp)
      complex*16, intent(out) :: lambda_stat(n,n,nsp)
      integer :: ispin
      complex*16, allocatable :: vx(:,:)
      do ispin=1,nsp
        allocate(vx(n,n))
        vx=(0.d0,0.d0)
        if(irel<2) call vertex_vk_0_r(vloc,k_pw(1,1,ispin),vx,
     &                                fif(1,1,1,1,ispin),n,np,
     &                                n_pbmtm_red)
        if(irel==2) call vertex_vk_0_c(vloc,k_pw,vx,fif,n,np,
     &                                 n_pbmtm_red)
        lambda_stat(:,:,ispin)=vx
        deallocate(vx)
      enddo  !! over ispin
      end
