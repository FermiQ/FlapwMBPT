      subroutine sigc_ii_mem(ind_tau,sig_tmp)
c     ------ (G_tau)*(W_tau) contribution to Sig_C -----------------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau
      complex*16, intent(inout) :: sig_tmp(nbndf,nbndf,2,
     &                                     ndim3_k(me_k+1),nspin)
      integer :: ind,i,i0,j,j0,ind_k,k,ig,ir,k0,irl,jrl,ispin,it,ista,
     &           ind_ista,itld,iq,iq0,ii,jj,nbas,iba,krl,jba,nb,ind_i,
     &           nbnd,k1,j00,ivs,i00,nnn,npw,j01,ia,ib,ic,nst,n0,n1,n2
      real*8 :: pi2,v(3),phase0,phase,dd
      complex*16 :: cc,yy
      integer, allocatable :: n0_b(:),nd_b(:)
      real*8, allocatable :: wr(:,:),sig_r(:,:)
      complex*16, allocatable :: xx(:,:),ex(:),g_c(:,:,:,:,:),
     &                           w_c(:,:,:,:),wq(:,:),gk(:,:,:),
     &                           sk(:,:,:,:),g0(:,:),
     &                           sig(:,:,:,:,:),ff(:,:),sig1(:,:,:,:),
     &                           sig2(:,:),sig3(:,:),sig_c(:,:,:,:),
     &                           tmp(:,:),w0(:)
      pi2=pi+pi
      dd=sqrt(amega)/nr_full
      nnn=max(nbasmpw,nplw_gw*nrel)
      if(ncmpl==2) then
        allocate(w_c(nplw_gw,npnt,max_r_star,ndim4_nr_red(me_k_b+1)))
        call w_q_g_c1(ind_tau,w_c)
      else
        allocate(w_c(nplw_gw,npnt,max_r_inv,ndim4_nr_red(me_k_b+1)))
        call w_q_g_r1(ind_tau,w_c)
      endif
      if(ind_tau==1) call timel('** W_Q_G_R1 finished ***************')
      if(complex_ro.or.irel==2) then
        allocate(g_c(nnn,npnt,max_r_star,ndim4_nr_red(me_k_b+1),nrel))
      else
        allocate(g_c(nnn,npnt,max_r_inv,ndim4_nr_red(me_k_b+1),nrel))
      endif
c ------------------------------------------------------------------
      allocate(xx(nrel,nrel))
      allocate(ex(nr_full))
      allocate(nd_b(nproc_b))
      allocate(n0_b(nproc_b))
      call size_shift_par(nplw_gw,nproc_b,nd_b,n0_b)
      n0=n0_b(me_b+1)
      n1=n0+nd_b(me_b+1)
      do ispin=1,nspin
        do it=1,2
          if(complex_ro.or.irel==2) then
            call g_k_g_c1(it,ind_tau,g_c,ispin,nnn)
          else
            call g_k_g_r1(it,ind_tau,g_c,ispin,nnn)
          endif
c ------------------------------------------------------------------
          do ind_ista=1,ndim4_nr_red(me_k_b+1)
            ista=n4_mpi_nr_red(me_k_b+1)+ind_ista
            j=r_pnt_ind(ista)
            allocate(wq(nplw_gw,nqdiv))
            allocate(gk(nbasmpw,nrel,nqdiv))
            do iq=1,nqdiv
              ig=k_group(iq)
              iq0=i_kref(iq)
              npw=nplwgw(iq0)
              allocate(w0(npw))
              j0=ip_r(j,ig)
              j00=ip_r_star(1,ista,ig)
              if(complex_ro) then
                w0=w_c(1:npw,iq0,j00,ind_ista)
              else
                j01=ip_r_inv(j00,ista)
                if(j01>0) then
                  w0=w_c(1:npw,iq0,j01,ind_ista)
                else
                  i=ind_r_star(j00,ista)
                  phase=-pi2*dot_product(pnt(:,iq0),
     &                                   tshift_r(:,i,inv_num))
                  cc=dcmplx(cos(phase),sin(phase))
                  w0=cc*conjg(w_c(1:npw,iq0,-j01,ind_ista))
                endif
              endif
              v=shift(:,ig)+tshift_r(:,j,ig)
              phase0=pi2*dot_product(pnt(:,iq),v)
              itld=g_sym_0(iq,ig)
              do i=1,npw
                ii=indpw_gw(i,iq)
                ia=igbs(1,ii)+igbs(1,itld)
                ib=igbs(2,ii)+igbs(2,itld)
                ic=igbs(3,ii)+igbs(3,itld)
                jj=indplw(ia,ib,ic)
                jj=gbs_sym(jj,ig)
                jj=iplf_gk(jj,iq0)
                phase=phase0+pi2*dot_product(gbs(:,ii),shift(:,ig))
                cc=dcmplx(cos(phase),-sin(phase))
                wq(i,iq)=cc*w0(jj)
              enddo
              deallocate(w0)
              nbas=nbask(iq0)/nrel
              allocate(g0(nbask(iq0),nrel))
              if(complex_ro.or.irel==2) then
                do irl=1,nrel
                  g0(:,irl)=g_c(1:nbask(iq0),iq0,j00,ind_ista,irl)
                enddo
              else
                j01=ip_r_inv(j00,ista)
                if(j01>0) then
                  do irl=1,nrel
                    g0(:,irl)=g_c(1:nbask(iq0),iq0,j01,ind_ista,irl)
                  enddo
                else
                  i=ind_r_star(j00,ista)
                  phase=-pi2*dot_product(pnt(:,iq0),
     &                                   tshift_r(:,i,inv_num))
                  cc=dcmplx(cos(phase),sin(phase))
                  do irl=1,nrel
                    g0(:,irl)=cc*conjg(g_c(1:nbask(iq0),iq0,-j01,
     &                                     ind_ista,irl))
                  enddo
                endif
              endif
              do i=1,nbas
                ii=indgb(i,iq)
                ia=igbs(1,ii)+igbs(1,itld)
                ib=igbs(2,ii)+igbs(2,itld)
                ic=igbs(3,ii)+igbs(3,itld)
                jj=indplw(ia,ib,ic)
                jj=gbs_sym(jj,ig)
                jj=iplf_bk(jj,iq0)
                phase=phase0+pi2*dot_product(gbs(:,ii),shift(:,ig))
                cc=dcmplx(cos(phase),-sin(phase))
                if(irel/=2) then
                  gk(i,1,iq)=cc*g0(jj,1)
                else
                  xx=(0.d0,0.d0)
                  do irl=1,nrel
                    iba=(irl-1)*nbas+jj
                    do jrl=1,nrel
                      do krl=1,nrel
                        jba=(krl-1)*2+jrl
                        xx(irl,jrl)=xx(irl,jrl)+g0(iba,krl)*uj(jba,ig)
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
              deallocate(g0)
            enddo    !! over iq
c ------------------------------------------------------------------
            allocate(wr(nqdiv,nr_full))
            call w_rr_r_r1(wq,wr)
            deallocate(wq)
            allocate(sig_c(nqdiv,nr_full,nrel,nrel))
            if(complex_ro.or.irel==2) then
              call g_rr_r_c1(gk,sig_c)
            else
              allocate(sig_r(nqdiv,nr_full))
              call g_rr_r_r1(gk,sig_r)
            endif
            deallocate(gk)
            if(complex_ro.or.irel==2) then
              do jrl=1,nrel
                do irl=1,nrel
                  do i=1,nr_full
                    do ir=1,nqdiv
                      sig_c(ir,i,irl,jrl)=-sig_c(ir,i,irl,jrl)*wr(ir,i)
                    enddo
                  enddo   !! over i
                enddo  !! over irl
              enddo  !! over jrl
            else
              do i=1,nr_full
                do ir=1,nqdiv
                  sig_r(ir,i)=-sig_r(ir,i)*wr(ir,i)
                enddo
              enddo   !! over i
              sig_c(:,:,1,1)=sig_r
              deallocate(sig_r)
            endif
            deallocate(wr)
            call fft3(ndiv(1),ndiv(2),ndiv(3),nr_full*nrel**2,sig_c,-1)
c ------------------------------------------------------------------
            allocate(ff(nr_full,nrel))
            allocate(sk(nplw_gw,nrel,nrel,nqdiv))
            do k=1,nqdiv
              call zone1_number(pnt(1,k),rb0,ndiv,k1)
              k0=i_kref(k)
              nbnd=n_bnd(k0,ispin)
              nb=nbask(k0)
              nbas=nb/nrel
c ----- Phase factors ----------------------------------------------
              do i=1,nr_full
                v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                                   +r_pnt(3,i)*r0b0(:,3)
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
c ------------------------------------------------------------------
            deallocate(ff,sig_c)
            if(complex_ro.or.irel==2) then
              nst=num_r_star(ista)
            else
              nst=num_r_inv(ista)
            endif
            do k=1,npnt
              npw=nplwgw(k)
              do ivs=1,nst
                if(complex_ro.or.irel==2) then
                  ir=ind_r_star(ivs,ista)
                else
                  ir=ind_r_inv(ivs,ista)
                  ir=i_r_inv(ir)
                endif
                ig=r_group(ir)
                k0=k_sym_0(k,ig)
                v=shift(:,ig)+tshift_r(:,ir,ig)
                phase0=pi2*dot_product(pnt(:,k),v)
                itld=g_sym_0(k,ig)
                do i=1,npw
                  ii=indpw_gw(i,k)
                  ia=igbs(1,ii)+igbs(1,itld)
                  ib=igbs(2,ii)+igbs(2,itld)
                  ic=igbs(3,ii)+igbs(3,itld)
                  jj=indplw(ia,ib,ic)
                  jj=gbs_sym(jj,ig)
                  jj=iplf_gk(jj,k0)
                  phase=phase0+pi2*dot_product(gbs(:,ii),shift(:,ig))
                  cc=dcmplx(cos(phase),-sin(phase))
                  if(irel/=2) then
                    g_c(i,k,ivs,ind_ista,1)=cc*sk(jj,1,1,k0)
                  else
                    xx=(0.d0,0.d0)
                    do jrl=1,nrel
                      do krl=1,nrel
                        jba=(krl-1)*2+jrl
                        do irl=1,nrel
                          xx(irl,jrl)=xx(irl,jrl)+sk(jj,irl,krl,k0)
     &                                           *uj(jba,ig)
                        enddo
                      enddo
                    enddo
                    do irl=1,nrel
                      i00=(irl-1)*nplw_gw+i
                      do jrl=1,nrel
                        yy=(0.d0,0.d0)
                        do krl=1,nrel
                          jba=(krl-1)*2+irl
                          yy=yy+conjg(uj(jba,ig))*xx(krl,jrl)
                        enddo
                        g_c(i00,k,ivs,ind_ista,jrl)=cc*yy
                      enddo
                    enddo
                  endif
                enddo   !! over i
              enddo   !! over ivs
            enddo   !! over k
            deallocate(sk)
          enddo   !! over ind_ista
c ------------------------------------------------------------------
          if(complex_ro.or.irel==2) then
            allocate(sig(nr_full,nd_b(me_b+1),ndim3_k(me_k+1),nrel,
     &                   nrel))
            allocate(tmp(nr_full,npnt))
          else
            allocate(sig(nr_full_inv,nd_b(me_b+1),ndim3_k(me_k+1),nrel,
     &                   nrel))
            allocate(tmp(nr_full_inv,npnt))
          endif
          sig=(0.d0,0.d0)
          do jrl=1,nrel
            do irl=1,nrel
              do i=1,nplw_gw
                ii=(irl-1)*nplw_gw+i
                tmp=(0.d0,0.d0)
                do k=1,npnt
                  do ind_ista=1,ndim4_nr_red(me_k_b+1)
                    ista=n4_mpi_nr_red(me_k_b+1)+ind_ista
                    if(complex_ro.or.irel==2) then
                      do ivs=1,num_r_star(ista)
                        ir=ind_r_star(ivs,ista)
                        tmp(ir,k)=g_c(ii,k,ivs,ind_ista,jrl)
                      enddo
                    else
                      do ivs=1,num_r_inv(ista)
                        ir=ind_r_inv(ivs,ista)
                        tmp(ir,k)=g_c(ii,k,ivs,ind_ista,jrl)
                      enddo
                    endif
                  enddo
                enddo
                if(nproc_k*nproc_b/=1) then
                  if(complex_ro.or.irel==2) then
                    call dgop(tmp,2*nr_full*npnt,'  +',comm_k_b)
                  else
                    call dgop(tmp,2*nr_full_inv*npnt,'  +',comm_k_b)
                  endif
                endif
                if(i<=n0) cycle
                if(i>n1) cycle
                do ind_k=1,ndim3_k(me_k+1)
                  k=n3_mpi_k(me_k+1)+ind_k
                  sig(:,i-n0,ind_k,irl,jrl)=tmp(:,k)
                enddo
              enddo
            enddo
          enddo
          deallocate(tmp)
c ------------------------------------------------------------------
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            npw=nplwgw(k)
            nbas=nbask(k)/nrel
            nb=n_bnd(k,ispin)
            n2=min(nd_b(me_b+1),npw-n0)
            allocate(sig3(nb,nb))
            sig3=(0.d0,0.d0)
            if(n2<=0) goto 1
c ----- Phase factors ------------------
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                                 +r_pnt(3,i)*r0b0(:,3)
              phase=pi2*dot_product(pnt(:,k),v)
              ex(i)=dcmplx(cos(phase),sin(phase))
            enddo
            allocate(sig1(n2,npw,nrel,nrel))
            allocate(ff(nr_full,n2))
            if(complex_ro.or.irel==2) then
              do jrl=1,nrel
                do irl=1,nrel
                  do ind_i=1,n2
                    do ir=1,nr_full
                      ff(ir,ind_i)=ex(ir)*sig(ir,ind_i,ind_k,irl,jrl)
                    enddo
                  enddo
                  call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n2,ff,1)
                  do i=1,npw
                    i0=indpw_gw(i,k)
                    call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
                    do ind_i=1,n2
                      sig1(ind_i,i,irl,jrl)=ff(ind,ind_i)
                    enddo
                  enddo
                enddo
              enddo   !! over jrl
            else
              do ir=1,nr_full
                i=nr_full_ind(ir)
                if(i>0) then
                  cc=ex(ir)
                  do ind_i=1,n2
                    ff(ir,ind_i)=cc*sig(i,ind_i,ind_k,1,1)
                  enddo
                else
                  phase=-pi2*dot_product(pnt(:,k),
     &                                   tshift_r(:,ir,inv_num))
                  cc=dcmplx(cos(phase),sin(phase))*ex(ir)
                  do ind_i=1,n2
                    ff(ir,ind_i)=cc*conjg(sig(-i,ind_i,ind_k,1,1))
                  enddo
                endif
              enddo
              call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n2,ff,1)
              do i=1,npw
                i0=indpw_gw(i,k)
                call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
                do ind_i=1,n2
                  sig1(ind_i,i,1,1)=ff(ind,ind_i)
                enddo
              enddo
            endif
            deallocate(ff)
            sig1=sig1*dd
c ------------------------------------------------------------------
            allocate(sig2(n2,nb))
            do irl=1,nrel
              sig2=(0.d0,0.d0)
              do jrl=1,nrel
                call zgemm('n','n',n2,nb,npw,(1.d0,0.d0),
     &                     sig1(1,1,irl,jrl),n2,
     &                     sa(1,1,ind_k,jrl,ispin),nplw_gw,(1.d0,0.d0),
     &                     sig2,n2)
              enddo
              call zgemm('c','n',nb,nb,n2,(1.d0,0.d0),
     &                   sa(n0+1,1,ind_k,irl,ispin),nplw_gw,sig2,n2,
     &                   (1.d0,0.d0),sig3,nb)
            enddo
            deallocate(sig1,sig2)
1           continue
            if(nproc_b/=1) call dgop(sig3,2*nb*nb,'  +',comm_b)
            sig_tmp(1:nb,1:nb,it,ind_k,ispin)=
     &        sig_tmp(1:nb,1:nb,it,ind_k,ispin)+sig3
            deallocate(sig3)
          enddo  !! over ind_k
          deallocate(sig)
        enddo   !! over it
      enddo   !! over ispin
      deallocate(w_c,g_c,xx,ex,n0_b,nd_b)
      if(ind_tau==1) call timel('** SIGC_II_MEM finished ************')
      end
