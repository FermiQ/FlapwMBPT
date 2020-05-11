      subroutine k0_bnd_stat_tau(ispin,ind_tau,iq,k0_tau,g_tau,ps3,s2)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,s2,ind_tau,ispin
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c)
      complex*16, intent(in) :: ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,
     &                              nqdiv_c)
	  complex*16, intent(out) :: k0_tau(nbndf_bnd,nbndf_bnd,nqdiv_c,2)
      integer :: k,kq,n,k0,nk,kq0,nkq,it,kf0,kfq0,ka
      real*8 :: v(3)
      complex*16, allocatable :: tmp(:,:),gt(:,:,:,:)
      n=nbndf_bnd
      allocate(tmp(n,n))
c ------------------------------------------------------------------
      k0_tau=(0.d0,0.d0)
c ------- Unpacking G_TAU ------------------------------------------
      allocate(gt(n,n,2,npnt_c))
      gt=(0.d0,0.d0)
	  do k=1,npnt_c
        ka=k_a_from_c(k)
	    call ferm_unpack_tau(gt(1,1,1,k),g_tau(1,1,1,ind_tau,k),
     &	                     n_low_bnd(ka,ispin),n,n,3)
      enddo
	  do k=1,nqdiv_c
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
	    v=pnt_c(:,k)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K-Q
	    kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kfq0,ispin)
	    do it=1,2
	      call zgemm('n','n',nk,nkq,nkq,(1.d0,0.d0),ps3(1,1,s2,k),n,
     &               gt(1,1,3-it,kq0),n,(0.d0,0.d0),tmp,n)
          call zgemm('n','n',nk,nkq,nk,(1.d0,0.d0),gt(1,1,it,k0),n,tmp,
     &               n,(0.d0,0.d0),k0_tau(1,1,k,it),n)
        enddo
      enddo  !! over k
      deallocate(tmp,gt)
      end
