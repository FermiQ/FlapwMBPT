      subroutine k_real_space_sf(iq,k_pw,k_rs,z_red,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq,ispin
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd,nqdiv_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(out) :: k_rs(nfun_red,nfun_red,
     &                                ndimc_kk(me_k+1))
      integer :: k0,nk,kq,kq0,nkq,ind_k,ind_r,ir,k,ka0,kaq0
      real*8 :: v(3)
      complex*16, allocatable :: t(:,:),t1(:,:,:)
      allocate(t(nbndf_bnd,nfun_red))
      allocate(t1(nfun_red,nfun_red,nqdiv_c))
      t1=(0.d0,0.d0)
      do ind_k=1,ndimc_kk(me_k+1)
        k=nc_mpi_kk(me_k+1)+ind_k
        k0=i_kref_c(k)
        ka0=k_a_from_c(k0)
        nk=n_low_bnd(ka0,ispin)
        v=pnt_c(:,k)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)     !! for K-Q
        kq0=i_kref_c(kq)
        kaq0=k_a_from_c(kq0)
        nkq=n_low_bnd(kaq0,ispin)
        call zgemm('n','c',nk,nfun_red,nkq,(1.d0,0.d0),k_pw(1,1,k),
     &             nbndf_bnd,z_red(1,1,kq,ispin),nfun_red,(0.d0,0.d0),t,
     &             nbndf_bnd)
        call zgemm('n','n',nfun_red,nfun_red,nk,(1.d0,0.d0),
     &             z_red(1,1,k,ispin),nfun_red,t,nbndf_bnd,(0.d0,0.d0),
     &             t1(1,1,k),nfun_red)
      enddo
      if(nproc_k/=1) call dgop(t1,2*nfun_red**2*nqdiv_c,'  +',comm_k)
      call k_k_to_rr_all(t1,nfun_red**2)
      do ind_r=1,ndimc_kk(me_k+1)
        ir=nc_mpi_kk(me_k+1)+ind_r
        k_rs(:,:,ind_r)=t1(:,:,ir)
      enddo
      deallocate(t,t1)
      end
