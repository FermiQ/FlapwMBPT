      subroutine vertex_sf_to_bnd_mm(vx,lambda_tau,iq,z_red,ispin,ch)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      integer, intent(in) :: iq,ispin
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          vx(nfun_red,nfun_red,nqdiv_c)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        ndimc_kk(me_k+1))
      integer :: k0,nk,kq,kq0,k,nkq,kf0,kfq0,ind_k
      real*8 :: v(3)
      complex*16, allocatable :: tt(:,:)
      allocate(tt(nfun_red,nbndf_bnd))
      do ind_k=1,ndimc_kk(me_k+1)
        k=nc_mpi_kk(me_k+1)+ind_k
        k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
        nk=n_low_bnd(kf0,ispin)
        v=pnt_c(:,k)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)     !! for K-Q
        kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
        nkq=n_low_bnd(kfq0,ispin)
        if(ch=='ph') then
          call zgemm('n','n',nfun_red,nkq,nfun_red,(1.d0,0.d0),
     &               vx(1,1,k),nfun_red,z_red(1,1,kq,ispin),nfun_red,
     &               (0.d0,0.d0),tt,nfun_red)
          call zgemm('c','n',nk,nkq,nfun_red,(1.d0,0.d0),
     &               z_red(1,1,k,ispin),nfun_red,tt,nfun_red,
     &               (1.d0,0.d0),lambda_tau(1,1,ind_k),nbndf_bnd)
        else if(ch=='pp') then
          call zgemm('t','n',nfun_red,nkq,nfun_red,(1.d0,0.d0),
     &               vx(1,1,k),nfun_red,z_red(1,1,kq,ispin),nfun_red,
     &               (0.d0,0.d0),tt,nfun_red)
          call zgemm('c','n',nk,nkq,nfun_red,(1.d0,0.d0),
     &               z_red(1,1,k,ispin),nfun_red,tt,nfun_red,
     &               (1.d0,0.d0),lambda_tau(1,1,ind_k),nbndf_bnd)
        endif
      enddo   !! over k
      deallocate(tt)
      end
