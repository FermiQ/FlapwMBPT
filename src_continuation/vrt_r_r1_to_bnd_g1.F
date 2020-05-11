      subroutine vrt_r_r1_to_bnd_g1(kq,k,nkq,na,z_red,a_red,phf,vx_r_r1,
     &                              vx_bnd_g1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nkq,na,kq,k
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd),
     &                          a_red(nbasmpw_red,nbndf_bnd),
     &                          phf(nr_full_red,nqdiv_c),
     &                          vx_r_r1(nrs_red,nrs_red)
      complex*16, intent(out) :: vx_bnd_g1(nkq,na)
      integer :: ka0,kqa0,kqa,ka,irl,ind,ir,ii,i,i0,j,j0,nbas,nbas1
      complex*16 :: cc
      complex*16, allocatable :: tt(:,:),t(:),t1(:,:)
      ka=k_a_from_c(k)
      ka0=i_kref(ka)
      kqa=k_a_from_c(kq)
      kqa0=i_kref(kqa)
      nbas=nbask_red(ka0)/nrel
      nbas1=nbask_red(kqa0)/nrel
      allocate(t(nr_full_red))
      allocate(t1(nbask_red(kqa0),nrs_red))
      do j=1,nrs_red
        t1=(0.d0,0.d0)
        do irl=1,nrel
          ind=nfun_red+(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=conjg(phf(ir,kq))
            t(ir)=cc*vx_r_r1(ind+ir,j)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t,-1)
          ii=(irl-1)*nbas1
          do i=1,nbas1
            i0=indgb_red(i,kqa)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(ii+i,j)=t(ind)
          enddo
        enddo   !! over irl
      enddo
      t1=t1/nr_full_red*sqrt(amega)
      allocate(tt(nkq,nrs_red))
      call zgemm('c','n',nkq,nrs_red,nbask_red(kqa0),(1.d0,0.d0),
     &           a_red,nbasmpw_red,t1,nbask_red(kqa0),(0.d0,0.d0),tt,
     &           nkq)
      deallocate(t1)
      call zgemm('c','n',nkq,nrs_red,nfun_red,(1.d0,0.d0),z_red,
     &           nfun_red,vx_r_r1,nrs_red,(1.d0,0.d0),tt,nkq)
      vx_bnd_g1(:,1:nfun_red)=tt(:,1:nfun_red)
      allocate(t1(nkq,nbas))
      do i=1,nkq
        t1=(0.d0,0.d0)
        do irl=1,nrel
          ind=nfun_red+(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=phf(ir,k)
            t(ir)=cc*tt(i,ind+ir)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t,1)
          ii=(irl-1)*nbas
          do j=1,nbas
            j0=indgb_red(j,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(i,ii+j)=t(ind)
          enddo
        enddo   !! over irl
      enddo
      t1=t1/nr_full_red*sqrt(amega)
      vx_bnd_g1(:,nfun_red+1:na)=t1
      deallocate(tt,t,t1)
      end
