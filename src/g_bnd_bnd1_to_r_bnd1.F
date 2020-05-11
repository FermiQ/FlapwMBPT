      subroutine g_bnd_bnd1_to_r_bnd1(kq,nkq,z_red,a_red,phf,g_bnd_bnd1,
     &                                g_r_bnd1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nkq,kq
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd),
     &                          a_red(nbasmpw_red,nbndf_bnd),
     &                          phf(nr_full_red,nqdiv_c),
     &                          g_bnd_bnd1(nkq,nkq)
      complex*16, intent(out) :: g_r_bnd1(nrs_red,nkq)
      integer :: kqa0,kqa,irl,ind,ir,ii,i,i0,nbas1,j
      complex*16 :: cc
      complex*16, allocatable :: t(:),t1(:,:)
      kqa=k_a_from_c(kq)
      kqa0=i_kref(kqa)
      nbas1=nbask_red(kqa0)/nrel
      allocate(t1(nbask_red(kqa0),nkq))
      call zgemm('n','n',nfun_red,nkq,nkq,(1.d0,0.d0),z_red,nfun_red,
     &           g_bnd_bnd1,nkq,(0.d0,0.d0),g_r_bnd1,nrs_red)
      call zgemm('n','n',nbask_red(kqa0),nkq,nkq,(1.d0,0.d0),a_red,
     &           nbasmpw_red,g_bnd_bnd1,nkq,(0.d0,0.d0),t1,
     &           nbask_red(kqa0))
      allocate(t(nr_full_red))
      do j=1,nkq
        do irl=1,nrel
          ii=(irl-1)*nbas1
          t=(0.d0,0.d0)
          do i=1,nbas1
            i0=indgb_red(i,kqa)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t(ind)=t1(ii+i,j)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t,1)
          ind=nfun_red+(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=phf(ir,kq)/sqrt(amega)
            g_r_bnd1(ind+ir,j)=cc*t(ir)
          enddo
        enddo   !! over irl
      enddo
      deallocate(t,t1)
      end
