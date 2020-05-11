      subroutine sigx_ii_mem(ispin)
c     ------ (G_tau)*(W_tau) contribution to Sig_C -----------------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: ispin
      integer :: ind,i,i0,j,j0,ind_k,k,ig,ir,k0,irl,jrl,ista,ind_ista,
     &           itld,gbs_number,npr,iq,iq0,ii,jj,nbas,iba,krl,j00,
     &           jba,nb,ind_i,nbnd,ibas0,k1,khbr,ivs,i00,jbas0,ia,ib,ic,
     &           npw
      real*8 :: pi2,v(3),phase0,phase
      complex*16 :: cc,yy
      integer, allocatable :: md(:),nd(:)
      complex*16, allocatable :: xx(:,:),ex(:),g_c(:,:,:,:,:),
     &                           v_c(:,:,:,:),wq(:,:),gk(:,:,:),
     &                           sk(:,:,:,:),wr(:,:),sig(:,:,:,:),
     &                           ff(:,:),sig3(:,:),sig_c(:,:,:,:),
     &                           tmp(:,:),s(:,:),yyy(:,:,:),sig1(:,:),
     &                           sig2(:,:),ft(:,:,:)
      pi2=pi+pi
      khbr=0
      if(ubi==' hf'.and.hybrid_type>2) khbr=1
      if(ndima_nr_red(me+1)/=0) then
        allocate(v_c(nplw_gw,npnt,max_r_star,ndima_nr_red(me+1)))
      else
        allocate(v_c(1,1,1,1))
      endif
      call v_q_g_r1(v_c,khbr)
      if(ndima_nr_red(me+1)/=0) then
        allocate(g_c(nbasmpw,npnt,max_r_star,ndima_nr_red(me+1),nrel))
      else
        allocate(g_c(1,1,1,1,1))
      endif
      call gx_k_g_r1(g_c,ispin)
c ------------------------------------------------------------------
      allocate(xx(nrel,nrel))
      allocate(ex(nr_full))
c ------------------------------------------------------------------
      do ind_ista=1,ndima_nr_red(me+1)
        ista=na_mpi_nr_red(me+1)+ind_ista
        j=r_pnt_ind(ista)
        allocate(wq(nplw_gw,nqdiv))
        allocate(gk(nbasmpw,nrel,nqdiv))
        do iq=1,nqdiv
          ig=k_group(iq)
          iq0=i_kref(iq)
          j0=ip_r(j,ig)
          j00=ip_r_star(1,ista,ig)
          v=shift(:,ig)+tshift_r(:,j,ig)
          phase0=pi2*dot_product(pnt(:,iq),v)
          itld=g_sym_0(iq,ig)
          do i=1,nplwgw(iq0)
            ii=indpw_gw(i,iq)
            v=gbs(:,ii)+gbs(:,itld)
            jj=gbs_number(v)
            jj=gbs_sym(jj,ig)
            jj=iplf_gk(jj,iq0)
            phase=phase0+pi2*dot_product(gbs(:,ii),shift(:,ig))
            cc=dcmplx(cos(phase),-sin(phase))
            wq(i,iq)=cc*v_c(jj,iq0,j00,ind_ista)
          enddo
          nbas=nbask(iq0)/nrel
          do i=1,nbas
            ii=indgb(i,iq)
            v=gbs(:,ii)+gbs(:,itld)
            jj=gbs_number(v)
            jj=gbs_sym(jj,ig)
            jj=iplf_bk(jj,iq0)
            phase=phase0+pi2*dot_product(gbs(:,ii),shift(:,ig))
            cc=dcmplx(cos(phase),-sin(phase))
            if(irel/=2) then
              gk(i,1,iq)=cc*g_c(jj,iq0,j00,ind_ista,1)
            else
              xx=(0.d0,0.d0)
              do irl=1,nrel
                iba=(irl-1)*nbas+jj
                do jrl=1,nrel
                  do krl=1,nrel
                    jba=(krl-1)*2+jrl
                    xx(irl,jrl)=xx(irl,jrl)
     &                         +g_c(iba,iq0,j00,ind_ista,krl)*uj(jba,ig)
                  enddo
                enddo
              enddo
              do irl=1,nrel
                iba=(irl-1)*nbas+i
                do jrl=1,nrel
                  yy=(0.d0,0.d0)
                  do krl=1,nrel
                    jba=(krl-1)*2+irl
                    yy=yy+conjg(uj(jba,ig))*xx(krl,jrl)
                  enddo
                  gk(iba,jrl,iq)=cc*yy
                enddo
              enddo
            endif
          enddo
        enddo    !! over iq
c ------------------------------------------------------------------
        allocate(wr(nqdiv,nr_full))
        call w_rr_r_r1_x(wq,wr)
        deallocate(wq)
        allocate(sig_c(nqdiv,nr_full,nrel,nrel))
        call g_rr_r_r1_x(gk,sig_c)
        deallocate(gk)
        do jrl=1,nrel
          do irl=1,nrel
            do i=1,nr_full
              do ir=1,nqdiv
                sig_c(ir,i,irl,jrl)=-sig_c(ir,i,irl,jrl)*wr(ir,i)
              enddo
              call fft3(ndiv(1),ndiv(2),ndiv(3),1,sig_c(1,i,irl,jrl),-1)
            enddo   !! over i
          enddo  !! over irl
        enddo  !! over jrl
        deallocate(wr)
c ------------------------------------------------------------------
        allocate(ff(nr_full,nrel))
        allocate(sk(nplw_gw,nrel,nrel,nqdiv))
        do k=1,nqdiv
          call zone1_number(pnt(1,k),rb0,ndiv,k1)
          k0=i_kref(k)
          nbnd=n_bnd(k0,ispin)
          nb=nbask(k0)
          nbas=nb/nrel
c ----- Phase factors ------------------
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                               +r_pnt(3,i)*r0b0(:,3)
            phase=pi2*dot_product(pnt(:,k),v)
            ex(i)=dcmplx(cos(phase),sin(phase))
          enddo
          do irl=1,nrel
            do jrl=1,nrel
              do ir=1,nr_full
                ff(ir,jrl)=conjg(ex(ir))*sig_c(k1,ir,irl,jrl)
              enddo
            enddo
            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nrel,ff,-1)
            do i=1,nplwgw(k0)
              i0=indpw_gw(i,k)
              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
              do jrl=1,nrel
                sk(i,irl,jrl,k)=ff(ind,jrl)/nr_full*sqrt(amega)
              enddo
            enddo
          enddo   !! over irl
        enddo   !! over k
        deallocate(ff,sig_c)
        do k=1,npnt
          nbas=nbask(k)/nrel
          npw=nplwgw(k)
          allocate(s(nbas,npw))
          do j=1,npw
            j0=indpw_gw(j,k)
            do i=1,nbas
              i0=indgb(i,k)
              ia=igbs(1,j0)-igbs(1,i0)
              ib=igbs(2,j0)-igbs(2,i0)
              ic=igbs(3,j0)-igbs(3,i0)
              ind=indplw(ia,ib,ic)
              if(complex_ro) then
                s(i,j)=dcmplx(sovr(ind),sovi(ind))
              else
                s(i,j)=dcmplx(sovr(ind),0.d0)
              endif
            enddo
          enddo
          allocate(yyy(npw,nrel,nrel))
          allocate(ft(nbas,nrel,nrel))
          do ivs=1,num_r_star(ista)
            ir=ind_r_star(ivs,ista)
            ig=r_group(ir)
            k0=k_sym_0(k,ig)
            v=shift(:,ig)+tshift_r(:,ir,ig)
            phase0=pi2*dot_product(pnt(:,k),v)
            itld=g_sym_0(k,ig)
            yyy=(0.d0,0.d0)
            do i=1,npw
              ii=indpw_gw(i,k)
              v=gbs(:,ii)+gbs(:,itld)
              jj=gbs_number(v)
              jj=gbs_sym(jj,ig)
              jj=iplf_gk(jj,k0)
              phase=phase0+pi2*dot_product(gbs(:,ii),shift(:,ig))
              cc=dcmplx(cos(phase),-sin(phase))
              if(irel/=2) then
                yyy(i,1,1)=cc*sk(jj,1,1,k0)
              else
                xx=(0.d0,0.d0)
                do jrl=1,nrel
                  do krl=1,nrel
                    jba=(krl-1)*2+jrl
                    do irl=1,nrel
                      xx(irl,jrl)=xx(irl,jrl)+sk(jj,irl,krl,k0)
     &                                       *uj(jba,ig)
                    enddo
                  enddo
                enddo
                do irl=1,nrel
                  do jrl=1,nrel
                    yy=(0.d0,0.d0)
                    do krl=1,nrel
                      jba=(krl-1)*2+irl
                      yy=yy+conjg(uj(jba,ig))*xx(krl,jrl)
                    enddo
                    yyy(i,irl,jrl)=cc*yy
                  enddo
                enddo
              endif
            enddo   !! over i
            call zgemm('n','n',nbas,nrel*nrel,npw,(1.d0,0.d0),s,nbas,
     &                 yyy,npw,(0.d0,0.d0),ft,nbas)
            do i=1,nbas
              do irl=1,nrel
                i00=(irl-1)*nbas+i
                do jrl=1,nrel
                  g_c(i00,k,ivs,ind_ista,jrl)=ft(i,irl,jrl)
                enddo
              enddo
            enddo   !! over i
          enddo   !! over ivs
          deallocate(s,yyy,ft)
        enddo   !! over k
        deallocate(sk)
      enddo   !! over ind_ista
      if(ndima_nr_red(me+1)/=0) deallocate(v_c)
      allocate(sig(nr_full,nbasmpw,ndim3_k(me_k+1),nrel))
      allocate(tmp(nr_full,npnt))
      do jrl=1,nrel
        do i=1,nbasmpw
          tmp=(0.d0,0.d0)
          do k=1,npnt
            if(ndima_nr_red(me+1)/=0) then
              do ind_ista=1,ndima_nr_red(me+1)
                ista=na_mpi_nr_red(me+1)+ind_ista
                do ivs=1,num_r_star(ista)
                  ir=ind_r_star(ivs,ista)
                  tmp(ir,k)=g_c(i,k,ivs,ind_ista,jrl)
                enddo
              enddo
            endif
          enddo
#ifdef MPI
          if(nproc/=1) call dgop(tmp,2*nr_full*npnt,'  +',
     &                           MPI_COMM_WORLD)
#endif
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            sig(:,i,ind_k,jrl)=tmp(:,k)
          enddo
        enddo
      enddo
      deallocate(tmp)
      if(ndima_nr_red(me+1)/=0) deallocate(g_c)
c ------------------------------------------------------------------
      npr=nproc_t*nproc_b
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        nb=nbask(k)
        nbas=nb/nrel
        npw=nplwgw(k)
        allocate(s(nbas,npw))
        do j=1,npw
          j0=indpw_gw(j,k)
          do i=1,nbas
            i0=indgb(i,k)
            ia=igbs(1,j0)-igbs(1,i0)
            ib=igbs(2,j0)-igbs(2,i0)
            ic=igbs(3,j0)-igbs(3,i0)
            ind=indplw(ia,ib,ic)
            if(complex_ro) then
              s(i,j)=dcmplx(sovr(ind),sovi(ind))
            else
              s(i,j)=dcmplx(sovr(ind),0.d0)
            endif
          enddo
        enddo
        allocate(nd(npr))
        allocate(md(npr))
        call size_shift_par(nb,npr,nd,md)
        ibas0=md(me_t_b+1)
c ----- Phase factors ------------------
        do i=1,nr_full
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                             +r_pnt(3,i)*r0b0(:,3)
          phase=pi2*dot_product(pnt(:,k),v)
          ex(i)=dcmplx(cos(phase),sin(phase))
        enddo
        if(nd(me_t_b+1)/=0) then
          allocate(ff(nr_full,nd(me_t_b+1)))
          allocate(sig1(nd(me_t_b+1),npw))
          allocate(sig2(nd(me_t_b+1),nbas))
        endif
        allocate(sig3(nb,nb))
        sig3=(0.d0,0.d0)
        do jrl=1,nrel
          jbas0=(jrl-1)*nbas
          if(nd(me_t_b+1)/=0) then
            do ind_i=1,nd(me_t_b+1)
              do ir=1,nr_full
                ff(ir,ind_i)=ex(ir)*sig(ir,ibas0+ind_i,ind_k,jrl)
              enddo
            enddo
            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nd(me_t_b+1),ff,1)
            do j=1,npw
              j0=indpw_gw(j,k)
              call zone1_number(gbs(:,j0),rbas,nrdiv,ind)
              do ind_i=1,nd(me_t_b+1)
                sig1(ind_i,j)=ff(ind,ind_i)/nr_full*sqrt(amega)
              enddo
            enddo
            call zgemm('n','c',nd(me_t_b+1),nbas,npw,(1.d0,0.d0),sig1,
     &                 nd(me_t_b+1),s,nbas,(0.d0,0.d0),sig2,
     &                 nd(me_t_b+1))
            do j=1,nbas
              do ind_i=1,nd(me_t_b+1)
                sig3(ibas0+ind_i,jbas0+j)=sig2(ind_i,j)
              enddo
            enddo
          endif
        enddo   !! over jrl
        if(npr/=1) call dgop(sig3,2*nb*nb,'  +',comm_t_b)
        pgg(1:nb,1:nb,ind_k)=pgg(1:nb,1:nb,ind_k)+sig3
        if(nd(me_t_b+1)/=0) deallocate(ff,sig1,sig2)
        deallocate(sig3,nd,md,s)
      enddo   !! over ind_k
      deallocate(xx,ex,sig)
      end
