      subroutine k_real_space_im(iq,k,kcom,k_rs,a_red,phf,isp)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,isp,k
      complex*16, intent(in) :: kcom(nbndf_bnd,nfun_red),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,
     &                                nspin),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: k_rs(nrel*nr_full_red,nfun_red)
      integer :: k0,nk,kf0,nbas,ka,irl,ind,i,i0
      complex*16, allocatable :: t(:,:),t1(:,:)
	  k0=i_kref_c(k)
      kf0=k_a_from_c(k0)
	  nk=n_low_bnd(kf0,isp)
      nbas=nbask_red(kf0)/nrel
      ka=k_a_from_c(k)
      allocate(t(nbas,nfun_red))
      allocate(t1(nr_full_red,nfun_red))
      do irl=1,nrel
        ind=(irl-1)*nbas+1
        call zgemm('n','n',nbas,nfun_red,nk,(1.d0,0.d0),
     &             a_red(ind,1,k,isp),nbasmpw_red,kcom,nbndf_bnd,
     &             (0.d0,0.d0),t,nbas)
        t1=(0.d0,0.d0)
        do i=1,nbas
          i0=indgb_red(i,ka)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          t1(ind,:)=t(i,:)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nfun_red,t1,1)
        ind=(irl-1)*nr_full_red
        do i=1,nr_full_red
          k_rs(ind+i,:)=phf(i,k)*t1(i,:)
        enddo
      enddo  !! over irl
      deallocate(t,t1)
      k_rs=k_rs/sqrt(amega)
      end
