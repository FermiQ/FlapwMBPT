      subroutine self_energy_1_stat(ispin,ind_tau,iq,g_tau,s2,v2_q,
     &                              vrt2_tau,aw,fif,z_red,a_red,phf)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq,s2,ispin,ind_tau
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c),
     &                      fif(maxel_red**2,n_pbmtm_red,nsort)
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red),
     &                          vrt2_tau(nbndf_bnd,nrs_red,nqdiv_c,2),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
	  complex*16, intent(inout) :: aw(nra_red,nra_red,nqdiv_c,0:n_tau)
      integer :: k,kq,n,k0,nk,kq0,nkq,it,j,i,iatom,isort,nn,ind2,ka,
     &           mx2,indp,nnn,i_tau,ii,ind,i2i,i2,kf0,kfq0,iqa,na,ka0,
     &           ir,i0,irl,j0,nbas
      real*8 :: v(3)
      complex*16 :: cc
      complex*16, allocatable :: gt(:,:,:,:),v_mt(:),tmp(:),
     &                           v_int(:),dd_tau(:,:,:),tt(:,:),
     &                           grt(:,:),sigt(:,:,:,:),t0(:),gg(:,:)
      mx2=maxel_red**2
      n=nbndf_bnd
      iqa=k_a_from_c(iq)
      allocate(sigt(nrs_red,nrs_red,nqdiv_c,2))
      sigt=(0.d0,0.d0)
c ------- Unpacking G_TAU ------------------------------------------
      allocate(gt(n,n,2,npnt_c))
	  do k=1,npnt_c
        ka=k_a_from_c(k)
	    call ferm_unpack_tau(gt(1,1,1,k),g_tau(1,1,1,ind_tau,k),
     &	                     n_low_bnd(ka,ispin),n,n,3)
      enddo
c -------- Get V -------------------------------------------
      allocate(v_mt(n_lem2_red))
      allocate(v_int(nr_full_red))
      v_mt=(0.d0,0.d0)
      v_int=(0.d0,0.d0)
	  do iatom=1,natom
	    isort=is(iatom)
	    n=lfunm_red(isort)
	    nn=n_pbmt0_red(isort)
	    ind2=io_lem2_red(iatom)
	    indp=iopb_red(iatom)
	    call dgemm('n','t',2,n*n,nn,1.d0,v2_q(s2,indp),2*n_pbtot_red,
     &	           fif(1,1,isort),mx2,0.d0,v_mt(ind2),2)
      enddo
      v_int=(0.d0,0.d0)
      do i=1,nplwgw_red(iqa)
        i0=indpw_gw_red(i,iq)
        call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
        v_int(ind)=v2_q(s2,n_pbmt_red+i)
      enddo
      call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,v_int,-1)
      do ir=1,nr_full_red
        cc=conjg(phf(ir,iq))/sqrt(amega)
        v_int(ir)=cc*v_int(ir)
      enddo
c ---------------------------------------------------------------
      nnn=n_lem2_red+nr_full_red
      allocate(tmp(nrs_red))
      do kq=1,nqdiv_c  ! K+Q here
	    kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kfq0,ispin)
	    v=pnt_c(:,kq)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,k)
	    k=index_k1_c(k)     !! for K here
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
        allocate(dd_tau(nnn,nrs_red,nkq))
	    allocate(grt(nrs_red,nkq))
	    do it=1,2
	      do j=1,nkq
	        tmp=vrt2_tau(j,:,kq,3-it)
	        do i=1,n_lem2_red
	          dd_tau(i,:,j)=tmp*v_mt(i)
            enddo
	        do i=1,nr_full_red
	          dd_tau(n_lem2_red+i,:,j)=tmp*v_int(i)
            enddo
          enddo
c --- Transform G_bnd_bnd' to G_(tL)r_bnd' -------------------------
          allocate(gg(nkq,nkq))
          gg=gt(1:nkq,1:nkq,it,kq0)
          call g_bnd_bnd1_to_r_bnd1(kq,nkq,z_red(1,1,kq),a_red(1,1,kq),
     &                              phf,gg,grt)
          deallocate(gg)
          do ii=1,nkq   !! over bands
            do j=1,nrs_red
	          do iatom=1,natom
	            isort=is(iatom)
	            n=lfunm_red(isort)
	            ind=io_lem_red(iatom)-1
	            ind2=io_lem2_red(iatom)-1
	            i2i=ind2
	            do i=1,n
	              do i2=1,n
	                i2i=i2i+1
	                sigt(ind+i,j,k,it)=sigt(ind+i,j,k,it)
     &	                              +grt(ind+i2,ii)*dd_tau(i2i,j,ii)
                  enddo
                enddo
              enddo
              do i=1,nr_full_red
                i2=nfun_red+i
	            sigt(i2,j,k,it)=sigt(i2,j,k,it)+grt(nfun_red+i,ii)
     &                                   *dd_tau(n_lem2_red+i,j,ii)
              enddo
            enddo
          enddo
        enddo
        deallocate(dd_tau,grt)
      enddo  !! over kq
      deallocate(gt,v_mt,v_int,tmp)
      allocate(t0(nr_full_red))
      do k=1,nqdiv_c
        k0=i_kref_c(k)
        ka0=k_a_from_c(k0)
        n=nbask_red(ka0)
        nbas=n/nrel
        na=nfun_red+n
        allocate(tt(na,na))
        allocate(grt(na,nrs_red))
        i_tau=me_t*ndim3_tau+ind_tau-1
	    do it=1,2
	      if(it==2) i_tau=n_tau-i_tau
	      if(it==2.and.i_tau==n_tau/2) cycle
          grt(1:nfun_red,:)=sigt(1:nfun_red,:,k,it)
          if(n/=0) then
c --- Transform SIG_r_r' to SIG_G_r' -------------------------------
            do j=1,nrs_red
              do irl=1,nrel
                ind=nfun_red+(irl-1)*nr_full_red
                do ir=1,nr_full_red
                  cc=conjg(phf(ir,k))
                  t0(ir)=cc*sigt(ind+ir,j,k,it)
                enddo
                call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t0,
     &                    -1)
                ii=nfun_red+(irl-1)*nbas
                do i=1,nbas
                  i0=indgb_red(i,ka)
                  call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
                  grt(ii+i,j)=t0(ind)/nr_full_red*sqrt(amega)
                enddo
              enddo   !! over irl
            enddo
c ------------------------------------------------------------------
          endif
          tt(:,1:nfun_red)=grt(:,1:nfun_red)
          if(n/=0) then
c --- Transform SIG_G_r' to SIG_G_G' -------------------------------
            do i=1,na
              do irl=1,nrel
                ind=nfun_red+(irl-1)*nr_full_red
                do ir=1,nr_full_red
                  cc=phf(ir,k)
                  t0(ir)=cc*grt(i,ind+ir)
                enddo
                call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t0,
     &                    1)
                ii=nfun_red+(irl-1)*nbas
                do j=1,nbas
                  j0=indgb_red(j,ka)
                  call zone1_number(gbs(:,j0),rbas,nrdiv_red,ind)
                  tt(i,ii+j)=t0(ind)/nr_full_red*sqrt(amega)
                enddo
              enddo   !! over irl
            enddo
c ------------------------------------------------------------------
          endif
	      aw(1:na,1:na,k,i_tau)=aw(1:na,1:na,k,i_tau)-wgt_c(iq)*tt
        enddo
        deallocate(tt,grt)
      enddo
      deallocate(sigt,t0)
      end
