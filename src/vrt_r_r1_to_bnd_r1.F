      subroutine vrt_r_r1_to_bnd_r1(k,nk,z_red,a_red,phf,vx_r_r1,
     &                              vx_bnd_r1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nk,k
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd),
     &                          a_red(nbasmpw_red,nbndf_bnd),
     &                          phf(nr_full_red,nqdiv_c),
     &                          vx_r_r1(nrs_red,nrs_red)
      complex*16, intent(out) :: vx_bnd_r1(nk,nrs_red)
      integer :: ka0,ka,irl,ind,ir,ii,i,i0,nbas,j
      complex*16 :: cc
      complex*16, allocatable :: t(:),t1(:,:)
      ka=k_a_from_c(k)
      ka0=i_kref(ka)
      nbas=nbask_red(ka0)/nrel
      allocate(t(nr_full_red))
      allocate(t1(nbask_red(ka0),nrs_red))
      do j=1,nrs_red
        t1=(0.d0,0.d0)
        do irl=1,nrel
          ind=nfun_red+(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=conjg(phf(ir,k))
            t(ir)=cc*vx_r_r1(ind+ir,j)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t,-1)
          ii=(irl-1)*nbas
          do i=1,nbas
            i0=indgb_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(ii+i,j)=t(ind)
          enddo
        enddo   !! over irl
      enddo
      deallocate(t)
      t1=t1/nr_full_red*sqrt(amega)
      call zgemm('c','n',nk,nrs_red,nbask_red(ka0),(1.d0,0.d0),
     &           a_red,nbasmpw_red,t1,nbask_red(ka0),(0.d0,0.d0),
     &           vx_bnd_r1,nk)
      deallocate(t1)
      call zgemm('c','n',nk,nrs_red,nfun_red,(1.d0,0.d0),z_red,
     &           nfun_red,vx_r_r1,nrs_red,(1.d0,0.d0),vx_bnd_r1,nk)
      end
