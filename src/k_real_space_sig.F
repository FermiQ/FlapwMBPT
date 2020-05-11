      subroutine k_rs_sig_mm(k_pw,k_rs,z_red,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      complex*16, intent(in) :: k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c)
      complex*16, intent(out) :: k_rs(nfun_red,n_pbmt_red,nqdiv_c)
      integer :: ir,k0,nk,ind,j,k,kf0
      complex*16, allocatable :: kr(:,:)
      k_rs=(0.d0,0.d0)
      do k=1,nqdiv_c
        k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
        nk=n_low_bnd(kf0,ispin)
        call zgemm('n','n',nfun_red,n_pbmt_red,nk,(1.d0,0.d0),
     &             z_red(1,1,k),nfun_red,k_pw(1,1,k),nbndf_bnd,
     &             (0.d0,0.d0),k_rs(1,1,k),nfun_red)
      enddo  !! over k
      allocate(kr(nqdiv_c,nfun_red))
      do j=1,n_pbmt_red
        do k=1,nqdiv_c
          call zone1_number(pnt_c(1,k),rb0_c,ndiv_c,ind)
          kr(ind,:)=k_rs(:,j,k)
        enddo
        call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nfun_red,kr,1)
        do ir=1,nqdiv_c
          k_rs(:,j,ir)=kr(ir,:)/nqdiv_c
        enddo
      enddo
      deallocate(kr)
      end
