      subroutine k_real_space_mi_ii(iq,k,k_pw,kcom,a_red,phf,isp,jsp)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,k,isp,jsp
      complex*16, intent(in) :: k_pw(nbndf_bnd,nbndf_bnd),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,
     &                                nspin),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: kcom(nbndf_bnd,nrel*nr_full_red)
      integer :: k0,nk,kq,kq0,nkq,kf0,kfq0,ka,nbas0,nbas,irl,i,i0,ind,
     &           i1,kqa
      real*8 :: v(3)
      complex*16, allocatable :: tm(:,:),t1(:,:)
	  k0=i_kref_c(k)
      ka=k_a_from_c(k)
      kf0=k_a_from_c(k0)
	  nk=n_low_bnd(kf0,isp)
	  v=pnt_c(:,k)-pnt_c(:,iq)
	  call zone1_number(v,rb0_c,ndiv_c,kq)
	  kq=index_k1_c(kq)     !! for K-Q
      kqa=k_a_from_c(kq)
	  kq0=i_kref_c(kq)
      kfq0=k_a_from_c(kq0)
	  nkq=n_low_bnd(kfq0,jsp)
      nbas0=nbask_red(kfq0)
      nbas=nbas0/nrel      
      allocate(tm(nk,nbas0))
      call zgemm('n','c',nk,nbas0,nkq,(1.d0,0.d0),k_pw,nbndf_bnd,
     &           a_red(1,1,kq,jsp),nbasmpw_red,(0.d0,0.d0),tm,nk)
      tm=tm/sqrt(amega)
      allocate(t1(nr_full_red,nk))
      do irl=1,nrel
        t1=(0.d0,0.d0)
        do i=1,nbas
          i0=indgb_red(i,kqa)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          i1=(irl-1)*nbas+i
          t1(ind,:)=tm(:,i1)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nk,t1,-1)
        ind=(irl-1)*nr_full_red
        do i=1,nr_full_red
          kcom(1:nk,ind+i)=conjg(phf(i,kq))*t1(i,:)
        enddo
      enddo  !! over irl
      deallocate(tm,t1)
      end
