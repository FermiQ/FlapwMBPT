      subroutine k_real_space_vmt(iq,k_pw,k_rs,z_red,ispin)
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
      complex*16, intent(out) :: k_rs(n_lem2_red)
      integer :: k0,nk,kq,kq0,nkq,k,ka0,kaq0,i,iatom,ij,ind,ind2,isort,
     &           j,n,nn
      real*8 :: v(3)
      complex*16, allocatable :: t(:,:),t1(:,:)
      allocate(t(nbndf_bnd,nfun_red))
      allocate(t1(maxel_red,maxel_red))
      k_rs=(0.d0,0.d0)
      do k=1,nqdiv_c
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
        do iatom=1,natom
          isort=is(iatom)
          n=lfunm_red(isort)
          nn=n*n
          ind=io_lem_red(iatom)
          ind2=io_lem2_red(iatom)-1
          call zgemm('n','n',n,n,nk,(1.d0,0.d0),z_red(ind,1,k,ispin),
     &               nfun_red,t(1,ind),nbndf_bnd,(0.d0,0.d0),t1,
     &               maxel_red)
          ij=ind2
          do j=1,n
            do i=1,n
              ij=ij+1
              k_rs(ij)=k_rs(ij)+t1(i,j)
            enddo
          enddo
        enddo
      enddo
      k_rs=k_rs/nqdiv_c
      deallocate(t,t1)
      end
