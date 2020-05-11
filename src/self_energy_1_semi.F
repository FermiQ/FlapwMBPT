      subroutine self_energy_1_semi(ispin,iq,g_tau,s2,vrt_nu,aw,fif,
     &                              z_red,a_red,phf,w2_q)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq,s2,ispin
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c),
     &                      fif(maxel_red**2,n_pbmtm_red,nsort)
      complex*16, intent(in) :: vrt_nu(nbndf_bnd,nrs_red,nqdiv_c,
     &                                 ndim3_nu),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          w2_q(n_pbtot_red,n_pbtot_red,ndim3_nu),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(inout) :: aw(nra_red,nra_red,nqdiv_c,0:n_tau)
      integer :: k,kq,n,k0,nk,kq0,nkq,it,j,i,ii,ind_tau,ind_nu,mx2,
     &           iatom,isort,nn,ind2,indp,nnn,ind,i2i,i2,i_tau,
     &           kf0,kfq0,ka,iqa,ka0,na,ir,i0,irl,j0,nbas
      real*8 :: v(3)
      complex*16 :: cc
      complex*16, allocatable :: gt(:,:,:,:,:),gg(:,:),t0(:),
     &                           dd_nu(:,:,:),
     &                           dd_tau(:,:,:,:,:),w_mt(:,:),
     &                           w_int(:,:),grt(:,:),sigt(:,:),
     &                           grt1(:,:),tt1(:,:)
      mx2=maxel_red**2
      n=nbndf_bnd
      iqa=k_a_from_c(iq)
      allocate(sigt(nrs_red,nrs_red))
c ------- Unpacking G_TAU ------------------------------------------
      allocate(gt(n,n,2,ndim3_tau,npnt_c))
      do k=1,npnt_c
        ka=k_a_from_c(k)
        do ind_tau=1,ndim3_tau
          call ferm_unpack_tau(gt(1,1,1,ind_tau,k),
     &	                       g_tau(1,1,1,ind_tau,k),
     &	                       n_low_bnd(ka,ispin),n,n,3)
        enddo
      enddo
c -------- Get W(q,nu) -------------------------------------------
      allocate(w_mt(n_lem2_red,ndim3_nu))
      allocate(w_int(nr_full_red,ndim3_nu))
      w_mt=(0.d0,0.d0)
      w_int=(0.d0,0.d0)
      do ind_nu=1,ndim3_nu
        do iatom=1,natom
          isort=is(iatom)
          n=lfunm_red(isort)
          nn=n_pbmt0_red(isort)
          ind2=io_lem2_red(iatom)
          indp=iopb_red(iatom)
          call dgemm('n','t',2,n*n,nn,1.d0,w2_q(s2,indp,ind_nu),
     &	           2*n_pbtot_red,fif(1,1,isort),mx2,0.d0,
     &               w_mt(ind2,ind_nu),2)
        enddo
        do i=1,nplwgw_red(iqa)
          i0=indpw_gw_red(i,iq)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          w_int(ind,ind_nu)=w2_q(s2,n_pbmt_red+i,ind_nu)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,
     &            w_int(1,ind_nu),-1)
        do ir=1,nr_full_red
          cc=conjg(phf(ir,iq))/sqrt(amega)
          w_int(ir,ind_nu)=cc*w_int(ir,ind_nu)
        enddo
      enddo
c ---------------------------------------------------------------
      nnn=n_lem2_red+nr_full_red
      allocate(t0(nr_full_red))
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
        ka0=k_a_from_c(k0)
        nn=nbask_red(ka0)
        nbas=nn/nrel
        na=nfun_red+nn
        allocate(dd_tau(nnn,nrs_red,2,ndim3_tau,nkq))
        allocate(dd_nu(nnn,nrs_red,ndim3_nu))
        do j=1,nkq
          do ind_nu=1,ndim3_nu
            do i=1,n_lem2_red
              dd_nu(i,:,ind_nu)=vrt_nu(j,:,kq,ind_nu)
     &                         *w_mt(i,ind_nu)
            enddo
            do i=1,nr_full_red
              dd_nu(n_lem2_red+i,:,ind_nu)=
     &	      vrt_nu(j,:,kq,ind_nu)*w_int(i,ind_nu)
            enddo
          enddo
          call nu_to_tau_half(dd_nu,dd_tau(1,1,1,1,j),nnn*nrs_red)
        enddo
        deallocate(dd_nu)
c ----------------------------------------------------------------
        allocate(grt(nrs_red,nkq))
        allocate(tt1(na,na))
        allocate(grt1(na,nrs_red))
        do ind_tau=1,ndim3_tau
          i_tau=me_t*ndim3_tau+ind_tau-1
          do it=1,2
            if(it==2) i_tau=n_tau-i_tau
            if(it==2.and.i_tau==n_tau/2) cycle
c --- Transform G_bnd_bnd' to G_(tL)r_bnd' -------------------------
            allocate(gg(nkq,nkq))
            gg=gt(1:nkq,1:nkq,it,ind_tau,kq0)
            call g_bnd_bnd1_to_r_bnd1(kq,nkq,z_red(1,1,kq),
     &                                a_red(1,1,kq),phf,gg,grt)
            deallocate(gg)
            sigt=(0.d0,0.d0)
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
                      sigt(ind+i,j)=sigt(ind+i,j)+grt(ind+i2,ii)
     &                             *dd_tau(i2i,j,3-it,ind_tau,ii)
                    enddo
                  enddo
                enddo
                do i=1,nr_full_red
                  i2=nfun_red+i
                  sigt(i2,j)=sigt(i2,j)+grt(i2,ii)
     &                      *dd_tau(n_lem2_red+i,j,3-it,ind_tau,ii)
                enddo
              enddo
            enddo
            grt1(1:nfun_red,:)=sigt(1:nfun_red,:)
            if(nn/=0) then
c --- Transform SIG_r_r' to SIG_G_r' -------------------------------
              do j=1,nrs_red
                do irl=1,nrel
                  ind=nfun_red+(irl-1)*nr_full_red
                  do ir=1,nr_full_red
                    cc=conjg(phf(ir,k))
                    t0(ir)=cc*sigt(ind+ir,j)
                enddo
                call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t0,
     &                    -1)
                ii=nfun_red+(irl-1)*nbas
                do i=1,nbas
                  i0=indgb_red(i,ka)
                  call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
                  grt1(ii+i,j)=t0(ind)/nr_full_red*sqrt(amega)
                enddo
              enddo   !! over irl
            enddo
c ------------------------------------------------------------------
            endif
            tt1(:,1:nfun_red)=grt1(:,1:nfun_red)
            if(nn/=0) then
c --- Transform SIG_G_r' to SIG_G_G' -------------------------------
              do i=1,na
                do irl=1,nrel
                  ind=nfun_red+(irl-1)*nr_full_red
                  do ir=1,nr_full_red
                    cc=phf(ir,k)
                    t0(ir)=cc*grt1(i,ind+ir)
                  enddo
                  call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t0,
     &                      1)
                  ii=nfun_red+(irl-1)*nbas
                  do j=1,nbas
                    j0=indgb_red(j,ka)
                    call zone1_number(gbs(:,j0),rbas,nrdiv_red,ind)
                    tt1(i,ii+j)=t0(ind)/nr_full_red*sqrt(amega)
                  enddo
                enddo   !! over irl
              enddo
c ------------------------------------------------------------------
            endif
	      aw(1:na,1:na,k,i_tau)=aw(1:na,1:na,k,i_tau)-wgt_c(iq)*tt1
          enddo
        enddo
        deallocate(dd_tau,grt,grt1,tt1)
      enddo  !! over kq
      deallocate(gt,w_mt,w_int,sigt,t0)
      end
