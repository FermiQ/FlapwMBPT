      subroutine k_rs_s_mi(krs0,k_rs,z_red,ispin)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      complex*16, intent(in) :: krs0(nbndf_bnd,nr_full_red,nqdiv_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c)
      complex*16, intent(out) :: k_rs(nfun_red,nr_full_red,nqdiv_c)
      integer :: ir,k0,nk,ind,j,k,kf0
      complex*16, allocatable :: kr(:,:),t(:,:)
      allocate(t(nr_full_red,nbndf_bnd))
      k_rs=(0.d0,0.d0)
      do k=1,nqdiv_c
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
	    call zgemm('n','n',nfun_red,nr_full_red,nk,(1.d0,0.d0),
     &             z_red(1,1,k),nfun_red,krs0(1,1,k),nbndf_bnd,
     &             (0.d0,0.d0),k_rs(1,1,k),nfun_red)
      enddo  !! over k
      deallocate(t)
      allocate(kr(nqdiv_c,nfun_red))
      do j=1,nr_full_red
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
      
      
      subroutine k_rs_s_mi_ii(k1,k_pw,krs0,phf,ispin)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: k1,ispin
      complex*16, intent(in) :: k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: krs0(nbndf_bnd,nr_full_red,nqdiv_c)
      integer :: k0,nk,kq,kq0,ind,k,kaq0,kf0,kaq,i,i0
      real*8 :: v(3)
      complex*16, allocatable :: t(:,:)
      allocate(t(nr_full_red,nbndf_bnd))
      krs0=(0.d0,0.d0)
      do k=1,nqdiv_c
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
	    v=pnt_c(:,k1)-pnt_c(:,k)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)  !! for K'-K
	    kq0=i_kref_c(kq)
        kaq=k_a_from_c(kq)
        kaq0=k_a_from_c(kq0)
        t=(0.d0,0.d0)
        do i=1,nplwgw_red(kaq0)
          i0=indpw_gw_red(i,kaq)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          t(ind,1:nk)=k_pw(1:nk,n_pbmt_red+i,k)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nk,t,1)
        do i=1,nr_full_red
          krs0(1:nk,i,k)=phf(i,kq)*t(i,1:nk)
        enddo
      enddo  !! over k
      deallocate(t)
      krs0=krs0/sqrt(amega)
      end 
      
      
      
      
      subroutine k_rs_s_im(k_pw,k_rs,a_red,phf,ispin)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      complex*16, intent(in) :: k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: k_rs(nrel*nr_full_red,n_pbmt_red,
     &                                nqdiv_c)
      integer :: ir,k0,nk,ind,j,k,kf0,ka,nbas,irl,i,i0
      complex*16, allocatable :: kr(:,:),t(:,:),t1(:,:)
      allocate(t(nbasmpw_red/nrel,n_pbmt_red))
      allocate(t1(nr_full_red,n_pbmt_red))
      k_rs=(0.d0,0.d0)
      do k=1,nqdiv_c
        ka=k_a_from_c(k)
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
        nbas=nbask_red(kf0)/nrel
        do irl=1,nrel
          ind=(irl-1)*nbas+1
          call zgemm('n','n',nbas,n_pbmt_red,nk,(1.d0,0.d0),
     &               a_red(ind,1,k),nbasmpw_red,k_pw(1,1,k),nbndf_bnd,
     &               (0.d0,0.d0),t,nbasmpw_red/nrel)
          t1=(0.d0,0.d0)
          do i=1,nbas
            i0=indgb_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(ind,:)=t(i,:)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),n_pbmt_red,
     &              t1,1)
          ind=(irl-1)*nr_full_red
          do i=1,nr_full_red
            k_rs(ind+i,:,k)=phf(i,k)*t1(i,:)
          enddo
        enddo  !! over irl
      enddo  !! over k
      k_rs=k_rs/sqrt(amega)
      deallocate(t,t1)
      allocate(kr(nqdiv_c,nrel*nr_full_red))
      do j=1,n_pbmt_red
        do k=1,nqdiv_c
          call zone1_number(pnt_c(1,k),rb0_c,ndiv_c,ind)
          kr(ind,:)=k_rs(:,j,k)
        enddo
	    call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nrel*nr_full_red,kr,1)
        do ir=1,nqdiv_c
          k_rs(:,j,ir)=kr(ir,:)/nqdiv_c
        enddo
      enddo
      deallocate(kr)
      end
      
      
      
      
      
      
      
      subroutine k_rs_s_ii(krs0,k_rs,a_red,phf,ispin)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      complex*16, intent(in) :: krs0(nbndf_bnd,nr_full_red,nqdiv_c),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: k_rs(nrel*nr_full_red,nr_full_red,
     &                                nqdiv_c)
      integer :: ir,k0,nk,ind,j,k,kf0,ka,irl,ii,nbas,i,i0
      complex*16, allocatable :: tm(:,:),kr(:,:),t(:,:),t1(:,:)
      allocate(tm(nbasmpw_red,nr_full_red))
      allocate(t(nr_full_red,nbndf_bnd))
      allocate(t1(nr_full_red,nr_full_red))
      k_rs=(0.d0,0.d0)
      do k=1,nqdiv_c
        ka=k_a_from_c(k)
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
        nbas=nbask_red(kf0)/nrel
	    nk=n_low_bnd(kf0,ispin)
        call zgemm('n','n',nbask_red(kf0),nr_full_red,nk,
     &             (1.d0,0.d0),a_red(1,1,k),nbasmpw_red,
     &             krs0(1,1,k),nbndf_bnd,(0.d0,0.d0),tm,nbasmpw_red)
        do irl=1,nrel
          ii=(irl-1)*nbas
          t1=(0.d0,0.d0)
          do i=1,nbas
            i0=indgb_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(ind,:)=tm(ii+i,:)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nr_full_red,
     &              t1,1)
          ind=(irl-1)*nr_full_red
          do i=1,nr_full_red
            k_rs(ind+i,:,k)=phf(i,k)*t1(i,:)
          enddo
        enddo   !! over irl
      enddo  !! over k
      deallocate(tm,t,t1)
      k_rs=k_rs/sqrt(amega)
      allocate(kr(nqdiv_c,nrel*nr_full_red))
      do j=1,nr_full_red
        do k=1,nqdiv_c
          call zone1_number(pnt_c(1,k),rb0_c,ndiv_c,ind)
          kr(ind,:)=k_rs(:,j,k)
        enddo
	    call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nrel*nr_full_red,kr,1)
        do ir=1,nqdiv_c
          k_rs(:,j,ir)=kr(ir,:)/nqdiv_c
        enddo
      enddo
      deallocate(kr)
      end
