      subroutine vertex_rs_to_bnd_mi(vx,lambda_tau,iq,z_red,ax,phf,isp,
     &                              jsp)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,isp,jsp
      complex*16, intent(in) :: phf(nr_full_red,nqdiv_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,
     &                             nspin),
     &                          vx(nfun_red,nrel*nr_full_red,nqdiv_c)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        nqdiv_c)
      integer :: k0,nk,kq,kq0,k,nkq,ka0,kqa0,kqa,irl,ind,ka,nbas,npr,
     &           ir,ii,i,i0
      real*8 :: v(3)
      complex*16 :: cc
      complex*16, allocatable :: tt(:,:),t(:,:),t1(:,:)
      npr=nrel*nplw_gw_red
      allocate(tt(nfun_red,nbndf_bnd))
      allocate(t(nr_full_red,nfun_red))
      allocate(t1(nfun_red,npr))
c -------- Transform Vrt(tau;nu) ---- > Vrt(omega;nu) ------------
	  do k=1,nqdiv_c
        ka=k_a_from_c(k)
	    k0=i_kref_c(k)
        ka0=k_a_from_c(k0)
	    nk=n_low_bnd(ka0,isp)
	    v=pnt_c(:,k)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K-Q
        kqa=k_a_from_c(kq)
	    kq0=i_kref_c(kq)
        kqa0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kqa0,jsp)
        nbas=nbask_red(kqa0)/nrel
        t1=(0.d0,0.d0)
        do irl=1,nrel
          ind=(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=phf(ir,kq)
            t(ir,:)=cc*vx(:,ind+ir,k)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nfun_red,t,1)
          ii=(irl-1)*nplwgw_red(kqa0)
          do i=1,nplwgw_red(kqa0)
            i0=indpw_gw_red(i,kqa)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(:,ii+i)=t(ind,:)
          enddo
        enddo   !! over irl
        t1=t1/nr_full_red*sqrt(amega)
        call zgemm('n','n',nfun_red,nkq,nrel*nplwgw_red(kqa0),
     &             (1.d0,0.d0),t1,nfun_red,ax(1,1,kq,jsp),npr,
     &             (0.d0,0.d0),tt,nfun_red)
        call zgemm('c','n',nk,nkq,nfun_red,(1.d0,0.d0),z_red(1,1,k,isp),
     &	           nfun_red,tt,nfun_red,(1.d0,0.d0),
     &             lambda_tau(1,1,k),nbndf_bnd)
      enddo   !! over k
      deallocate(tt,t,t1)
      end
