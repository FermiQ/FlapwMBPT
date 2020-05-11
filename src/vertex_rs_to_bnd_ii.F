      subroutine vertex_rs_to_bnd_ii(vx,lambda_tau,iq,ax,phf,isp,jsp)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,isp,jsp
      complex*16, intent(in) :: vx(nrel*nr_full_red,nrel*nr_full_red,
     &                             nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,
     &                             nspin),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        nqdiv_c)
      integer :: k0,nk,kq,kq0,k,nkq,ka0,nrr,nrp1,nrp,ka,kqa0,kqa,irl,
     &           ind,ir,ii,i,i0,np,np1,nrpp
      real*8 :: v(3)
      complex*16 :: cc
      complex*16, allocatable :: tt(:,:),t(:,:),t1(:,:),t2(:,:)
      nrr=nrel*nr_full_red
      nrpp=nrel*nplw_gw_red
      allocate(tt(nrpp,nbndf_bnd))
      allocate(t(nr_full_red,nrr))
      allocate(t1(nrr,nrpp))
      allocate(t2(nrpp,nrpp))
      nrr=nrel*nr_full_red
c -------- Transform Vrt(tau;nu) ---- > Vrt(omega;nu) ------------
	  do k=1,nqdiv_c
	    k0=i_kref_c(k)
        ka0=k_a_from_c(k0)
        ka=k_a_from_c(k)
	    nk=n_low_bnd(ka0,isp)
	    v=pnt_c(:,k)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K-Q
	    kq0=i_kref_c(kq)
        kqa0=k_a_from_c(kq0)
        kqa=k_a_from_c(kq)
	    nkq=n_low_bnd(kqa0,jsp)
        np=nplwgw_red(ka0)
        np1=nplwgw_red(kqa0)
        nrp=nrel*np
        nrp1=nrel*np1
        do irl=1,nrel
          ind=(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=phf(ir,kq)
            t(ir,:)=cc*vx(:,ind+ir,k)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nrr,t,1)
          ii=(irl-1)*np1
          do i=1,np1
            i0=indpw_gw_red(i,kqa)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(:,ii+i)=t(ind,:)
          enddo
        enddo   !! over irl
        t1=t1/nr_full_red*sqrt(amega)
        do irl=1,nrel
          ind=(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=conjg(phf(ir,k))
            t(ir,1:nrp1)=cc*t1(ind+ir,1:nrp1)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nrp1,t,-1)
          ii=(irl-1)*np
          do i=1,np
            i0=indpw_gw_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t2(ii+i,1:nrp1)=t(ind,1:nrp1)
          enddo
        enddo
        t2=t2/nr_full_red*sqrt(amega)
        call zgemm('n','n',nrp,nkq,nrp1,(1.d0,0.d0),t2,nrpp,
     &             ax(1,1,kq,jsp),nrpp,(0.d0,0.d0),tt,nrpp)
        call zgemm('c','n',nk,nkq,nrp,(1.d0,0.d0),ax(1,1,k,isp),nrpp,
     &             tt,nrpp,(1.d0,0.d0),lambda_tau(1,1,k),nbndf_bnd)
      enddo   !! over k
      deallocate(tt,t,t1,t2)
      end
