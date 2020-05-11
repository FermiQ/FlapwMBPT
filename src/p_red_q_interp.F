      subroutine p_red_q_interp
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ind_q,iq,k,ka,i,i0,j,j0,ind,ind_nu,k0,ka0,np,isy,ig
      real*8 :: v(3),phase,pi2
      complex*16 :: bb,cc,cc1
      real*8, allocatable :: tmp(:,:,:)
      complex*16, allocatable :: p(:,:),dd(:),mm(:,:,:),mi(:,:,:),
     &                           ii(:,:,:),pmi(:,:),f(:),pii(:,:)
      pi2=pi+pi
      if(allocated(p_vrt_nu_q)) deallocate(p_vrt_nu_q)
      allocate(p_vrt_nu_q(n_pbtot_red,n_pbtot_red,ndim3_nu,
     &                    ndim3_k(me_k+1)))
      p_vrt_nu_q=0.d0
      allocate(tmp(n_pbtot_red,n_pbtot_red,npnt_c))
      allocate(p(n_pbtot_red,n_pbtot_red))
      allocate(f(nr_full_red))
      allocate(dd(npnt))
      dd=(0.d0,0.d0)
      bb=(1.d0,0.d0)
      if(nqdiv_c/=1) then
        bb=bb+sum(bk_bz_red)
        do k=1,npnt
          dd(k)=sum(ckq_bz(:,k))
        enddo
      endif
      do ind_nu=1,ndim3_nu
        tmp=0.d0
        do ind_q=1,ndim_k_red(me_k+1)
          iq=n_mpi_k_red(me_k+1)+ind_q
          tmp(:,:,iq)=p_vertex_nu_q(:,:,ind_nu,ind_q)
        enddo
        if(nproc_k/=1) call dgop(tmp,n_pbtot_red**2*npnt_c,'  +',comm_k)
        allocate(mm(n_pbmt_red,n_pbmt_red,ndim3_k(me_k+1)))
        allocate(mi(n_pbmt_red,nr_full_red,ndim3_k(me_k+1)))
        allocate(ii(nr_full_red,nr_full_red,ndim3_k(me_k+1)))
        mm=(0.d0,0.d0)
        mi=(0.d0,0.d0)
        ii=(0.d0,0.d0)
        do k=1,nqdiv_c
          ka=k_a_from_c(k)
          k0=i_kref_c(k)
          ka0=k_a_from_c(k0)
          np=n_pbmt_red+nplwgw_red(ka0)
          call unpack_hermit(p,tmp(1,1,k0),np,n_pbtot_red,n_pbtot_red,
     &                       (0.d0,0.d0),(1.d0,0.d0))
          call sym_w_red(k,k0,p)
c ----------------------------- MM part ----------------------------
          do ind_q=1,ndim3_k(me_k+1)
            iq=n3_mpi_k(me_k+1)+ind_q
            if(k/=1) cc1=ckq_bz(k,iq)-bk_bz_red(k)
            if(k==1) cc1=bb-dd(iq)
            do j=1,n_pbmt_red
              do i=1,n_pbmt_red
                mm(i,j,ind_q)=mm(i,j,ind_q)+p(i,j)*cc1
              enddo
            enddo
          enddo
c ----------------------------- MI part ----------------------------
          allocate(pmi(n_pbmt_red,nr_full_red))
          do i=1,n_pbmt_red
            f=(0.d0,0.d0)
            do j=1,nplwgw_red(ka0)
              j0=indpw_gw_red(j,ka)
              call zone1_number(gbs(:,j0),rbas,nrdiv_red,ind)
              f(ind)=p(i,n_pbmt_red+j)
            enddo
            call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,f,-1)
            pmi(i,:)=f/sqrt(amega)
          enddo
          do j=1,nr_full_red
            v=r0b0_red(:,1)*r_pnt_red(1,j)
     &       +r0b0_red(:,2)*r_pnt_red(2,j)
     &       +r0b0_red(:,3)*r_pnt_red(3,j)
            phase=pi2*dot_product(pnt_c(:,k),v)
            cc=dcmplx(cos(phase),-sin(phase))
            pmi(:,j)=cc*pmi(:,j)
          enddo
          do ind_q=1,ndim3_k(me_k+1)
            iq=n3_mpi_k(me_k+1)+ind_q
            if(k/=1) cc1=ckq_bz(k,iq)-bk_bz_red(k)
            if(k==1) cc1=bb-dd(iq)
            do j=1,nr_full_red
              do i=1,n_pbmt_red
                mi(i,j,ind_q)=mi(i,j,ind_q)+pmi(i,j)*cc1
              enddo
            enddo
          enddo
          deallocate(pmi)
c ----------------------------- II part ----------------------------
          allocate(pmi(nplwgw_red(ka0),nr_full_red))
          do i=1,nplwgw_red(ka0)
            f=(0.d0,0.d0)
            do j=1,nplwgw_red(ka0)
              j0=indpw_gw_red(j,ka)
              call zone1_number(gbs(:,j0),rbas,nrdiv_red,ind)
              f(ind)=p(n_pbmt_red+i,n_pbmt_red+j)
            enddo
            call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,f,-1)
            pmi(i,:)=f/sqrt(amega)
          enddo
          do j=1,nr_full_red
            v=r0b0_red(:,1)*r_pnt_red(1,j)
     &       +r0b0_red(:,2)*r_pnt_red(2,j)
     &       +r0b0_red(:,3)*r_pnt_red(3,j)
            phase=pi2*dot_product(pnt_c(:,k),v)
            cc=dcmplx(cos(phase),-sin(phase))
            pmi(:,j)=cc*pmi(:,j)
          enddo
          allocate(pii(nr_full_red,nr_full_red))
          do j=1,nr_full_red
            f=(0.d0,0.d0)
            do i=1,nplwgw_red(ka0)
              i0=indpw_gw_red(i,ka)
              call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
              f(ind)=pmi(i,j)
            enddo
            call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,f,1)
            pii(:,j)=f/sqrt(amega)
          enddo
          deallocate(pmi)
          do i=1,nr_full_red
            v=r0b0_red(:,1)*r_pnt_red(1,i)
     &       +r0b0_red(:,2)*r_pnt_red(2,i)
     &       +r0b0_red(:,3)*r_pnt_red(3,i)
            phase=pi2*dot_product(pnt_c(:,k),v)
            cc=dcmplx(cos(phase),sin(phase))
            pii(i,:)=cc*pii(i,:)
          enddo
          do ind_q=1,ndim3_k(me_k+1)
            iq=n3_mpi_k(me_k+1)+ind_q
            if(k/=1) cc1=ckq_bz(k,iq)-bk_bz_red(k)
            if(k==1) cc1=bb-dd(iq)
            do j=1,nr_full_red
              do i=1,nr_full_red
                ii(i,j,ind_q)=ii(i,j,ind_q)+pii(i,j)*cc1
              enddo
            enddo
          enddo
          deallocate(pii)
        enddo   !! over k
c ---------------- Symmetrization ----------------------------------
c        allocate(pii(nr_full_red,nr_full_red))
c        do ind_q=1,ndim3_k(me_k+1)
c          iq=n3_mpi_k(me_k+1)+ind_q
c          pii=(0.d0,0.d0)
c          do isy=1,num_sym_k(iq)
c            ig=list_sym_k(isy,iq)
c            do j=1,nr_full_red
c              j0=ip_r_red(j,ig)
c              do i=1,nr_full_red
c                i0=ip_r_red(i,ig)
c                v=tshift_r_red(:,i,ig)-tshift_r_red(:,j,ig)
c                phase=pi2*dot_product(pnt(:,iq),v)
c                cc=dcmplx(cos(phase),sin(phase))
c                pii(i,j)=pii(i,j)+cc*ii(i0,j0,ind_q)
c              enddo
c            enddo
c          enddo
c          ii(:,:,ind_q)=pii/num_sym_k(iq)
c        enddo
c        deallocate(pii)
c ------------------------------------------------------------------
        do ind_q=1,ndim3_k(me_k+1)
          iq=n3_mpi_k(me_k+1)+ind_q
          np=n_pbmt_red+nplwgw_red(iq)
c ---------------------- MM part -----------------------------------
          do j=1,n_pbmt_red
            do i=1,n_pbmt_red
              p(i,j)=mm(i,j,ind_q)
            enddo
          enddo
c ---------------------- MI part -----------------------------------
          allocate(pmi(nr_full_red,n_pbmt_red))
          do j=1,nr_full_red
            v=r0b0_red(:,1)*r_pnt_red(1,j)
     &       +r0b0_red(:,2)*r_pnt_red(2,j)
     &       +r0b0_red(:,3)*r_pnt_red(3,j)
            phase=pi2*dot_product(pnt(:,iq),v)
            cc=dcmplx(cos(phase),sin(phase))
            do i=1,n_pbmt_red
              pmi(j,i)=mi(i,j,ind_q)*cc
            enddo
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),n_pbmt_red,
     &              pmi,1)
          do j=1,nplwgw_red(iq)
            j0=indpw_gw_red(j,iq)
            call zone1_number(gbs(:,j0),rbas,nrdiv_red,ind)
            do i=1,n_pbmt_red
              p(i,n_pbmt_red+j)=pmi(ind,i)*sqrt(amega)/nr_full_red
              p(n_pbmt_red+j,i)=conjg(p(i,n_pbmt_red+j))
            enddo
          enddo
          deallocate(pmi)
c ---------------------- II part -----------------------------------
          allocate(pmi(nr_full_red,nr_full_red))
          do j=1,nr_full_red
            v=r0b0_red(:,1)*r_pnt_red(1,j)
     &       +r0b0_red(:,2)*r_pnt_red(2,j)
     &       +r0b0_red(:,3)*r_pnt_red(3,j)
            phase=pi2*dot_product(pnt(:,iq),v)
            cc=dcmplx(cos(phase),sin(phase))
            do i=1,nr_full_red
              pmi(j,i)=ii(i,j,ind_q)*cc
            enddo
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nr_full_red,
     &              pmi,1)
          allocate(pii(nr_full_red,nplwgw_red(iq)))
          do j=1,nplwgw_red(iq)
            j0=indpw_gw_red(j,iq)
            call zone1_number(gbs(:,j0),rbas,nrdiv_red,ind)
            do i=1,nr_full_red
              pii(i,j)=pmi(ind,i)*sqrt(amega)/nr_full_red
            enddo
          enddo
          deallocate(pmi)
          allocate(pmi(nr_full_red,nplwgw_red(iq)))
          do i=1,nr_full_red
            v=r0b0_red(:,1)*r_pnt_red(1,i)
     &       +r0b0_red(:,2)*r_pnt_red(2,i)
     &       +r0b0_red(:,3)*r_pnt_red(3,i)
            phase=pi2*dot_product(pnt(:,iq),v)
            cc=dcmplx(cos(phase),-sin(phase))
            do j=1,nplwgw_red(iq)
              pmi(i,j)=pii(i,j)*cc
            enddo
          enddo
          deallocate(pii)
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),
     &              nplwgw_red(iq),pmi,-1)
          do i=1,nplwgw_red(iq)
            i0=indpw_gw_red(i,iq)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            do j=1,nplwgw_red(iq)
              p(n_pbmt_red+i,n_pbmt_red+j)=pmi(ind,j)
     &                                    *sqrt(amega)/nr_full_red
            enddo
          enddo
          deallocate(pmi)
          allocate(pii(n_pbtot_red,n_pbtot_red))
          allocate(pmi(n_pbtot_red,n_pbtot_red))
          pii=p
          p=(0.d0,0.d0)
          do isy=1,num_sym_k(iq)
            ig=list_sym_k(isy,iq)
            pmi=pii
            call sym_p_red(iq,ig,pmi)
            p=p+pmi
          enddo
          p=p/num_sym_k(iq)
          deallocate(pii,pmi)
          call pack_hermit(p,p_vrt_nu_q(1,1,ind_nu,ind_q),np,
     &                     n_pbtot_red,n_pbtot_red,0.d0,1.d0)
        enddo  !! over ind_q
        deallocate(mm,mi,ii)
      enddo   !! over ind_nu
      deallocate(p,tmp,f,dd)
      call timel('*** P_RED_Q_INTERP finished ********')
      end
