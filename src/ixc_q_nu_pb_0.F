      subroutine ixc_q_nu_pb_0(ch,ind_k,xi0_nu,xi_nu,n_ch0,
     &                         ind_ch0,n_ix0,ind_ix0,n_ch,resp_stat,
     &                         eps1,eps_tns,chi_mt,chi0_mt,chi_neig,mtp)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ch
      integer, intent(in) :: ind_k,n_ch0,ind_ch0(n_ch0),n_ix0,
     &	                     ind_ix0(n_ix0),n_ch
      real*8, intent(in) :: mtp(n_pbmt)
      real*8, intent(inout) :: resp_stat(npnt,n_ch)
      complex*16, intent(inout) :: xi0_nu(0:n_nu,n_line,n_ch),
     &	                           xi_nu(0:n_nu,n_line,n_ch),
     &	                           eps1(0:n_nu,n_line),
     &                             chi_mt(natom,ndim3_nu),
     &                             chi0_mt(natom,ndim3_nu),
     &                             chi_neig(max_sites,nsort,ndim3_nu),
     &                             eps_tns(3,3,0:n_nu)
      logical :: sing
      integer :: k,n_pbt,i,ind_nu,i_nu,info,iq_l,j,n,nn,nn2
      real*8 :: pi2a,q,q2,pi8,p00
      complex*16 :: cc,ets(3,3)
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: d(:),v_tmpr(:,:),ii_nu(:,:,:),
     &                       mi_nu(:,:,:,:),mm_nu(:,:,:)
      complex*16, allocatable :: chi0(:,:),work(:,:),thet_q(:,:,:),
     &                           v_tmp(:,:),chi(:,:),ixc0(:,:),tmp(:,:),
     &                           tmp1(:,:),a(:,:)
      k=n3_mpi_k(me_k+1)+ind_k
      iq_l=k_npnt_in_line(k)
      pi8=8.d0*pi
      pi2a=(pi+pi)/par
      q2=pi2a**2*dot_product(pnt(:,k),pnt(:,k))
      q=sqrt(q2)
      n_pbt=n_pbmt+nplwgw(k)
      nn2=ncmpl*nd_b_pbint(me_b+1)
      if(ch=='___00'.or.ch=='___ZZ') then
        n=1
        nn=n_pbt
      else if(ch=='00_ZZ') then
        n=2
        nn=2*n_pbt
      else if(ch=='___XY') then
        n=2
        nn=n_pbt
      endif
      allocate(ipiv(nn))
      allocate(work(nn,nn))
      allocate(chi0(nn,nn))
      allocate(tmp(n_pbt,n_pbt))
      allocate(tmp1(n_pbmt,n_pbmt))
      allocate(chi(nn,nn))
      allocate(ixc0(nn,nn))
      allocate(a(nn,nn))
      allocate(d(nn))
c -------------------------------------------------------------------
      allocate(thet_q(n_pbt,n_pbt,n_ixc_0))
      call theta_q_ss(k,n_pbt,thet_q)
      if(ch=='___00'.or.ch=='00_ZZ') then
        allocate(v_tmp(n_pbt,n_pbt))
        if(ncmpl==1) then
          allocate(v_tmpr(n_pbt,n_pbt))
          call v_coul_full_r(nplwgw(k),v_tmpr,ind_k,0)
          call pb_c_from_r(n_pbt,n_pbt,n_pbt,v_tmp,v_tmpr,pnt(1,k))
          deallocate(v_tmpr)
        else
          call v_coul_full(nplwgw(k),v_tmp,ind_k,0)
        endif
        if(k/=1) thet_q(:,:,1)=thet_q(:,:,1)+v_tmp
      endif
      if(n==1.or.ch=='___XY') then
        j=ind_ix0(1)
        ixc0(:,:)=thet_q(:,:,j)
        if(ch=='___XY') ixc0(:,:)=ixc0(:,:)/2 
      else
        j=ind_ix0(1)
        ixc0(1:nn/2,1:nn/2)=thet_q(:,:,j)
        j=ind_ix0(2)
        ixc0(nn/2+1:nn,nn/2+1:nn)=thet_q(:,:,j)
        j=ind_ix0(3)
        ixc0(nn/2+1:nn,1:nn/2)=thet_q(:,:,j)
        ixc0(1:nn/2,nn/2+1:nn)=thet_q(:,:,j)
      endif
      deallocate(thet_q)
c   ------- Transform P_FULL from TAU to NU representation for given Q -
      allocate(mm_nu(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,n))
      allocate(mi_nu(ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),ndim3_tn,n))
      allocate(ii_nu(nn2,ndim3_tn,n))
      if(ch=='___XY') then
        mm_nu(:,:,1)=resp_lda_tau_mm(:,:,ind_k,ind_ch0(1))
        mi_nu(:,:,:,1)=resp_lda_tau_mi(:,:,:,ind_k,ind_ch0(1))
        ii_nu(:,:,1)=resp_lda_tau_ii(:,:,ind_k,ind_ch0(1))
        mm_nu(:,:,2)=resp_lda_tau_mm(:,:,ind_k,ind_ch0(2))
        mi_nu(:,:,:,2)=resp_lda_tau_mi(:,:,:,ind_k,ind_ch0(2))
        ii_nu(:,:,2)=resp_lda_tau_ii(:,:,ind_k,ind_ch0(2))
        call tau_to_nu_ba3ut(mm_nu(1,1,1),ncmpl*ndim_pbmt(me_b+1),
     &                       ncmpl*ndim_pbmt(me_b+1),1)
        call tau_to_nu_ba3(mi_nu(1,1,1,1),ncmpl*nd_b_pbmt(me_b+1),
     &                     nplw_gw,ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),
     &                     1)
        call tau_to_nu_ba3ut(ii_nu(1,1,1),nn2,nn2,1)
        call tau_to_nu_ba3ut(mm_nu(1,1,2),ncmpl*ndim_pbmt(me_b+1),
     &                       ncmpl*ndim_pbmt(me_b+1),2)
        call tau_to_nu_ba3(mi_nu(1,1,1,2),ncmpl*nd_b_pbmt(me_b+1),
     &                     nplw_gw,ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),
     &                     2)
        call tau_to_nu_ba3ut(ii_nu(1,1,2),nn2,nn2,2)
      else
        do i=1,n
          mm_nu(:,:,i)=resp_lda_tau_mm(:,:,ind_k,ind_ch0(i))
          mi_nu(:,:,:,i)=resp_lda_tau_mi(:,:,:,ind_k,ind_ch0(i))
          ii_nu(:,:,i)=resp_lda_tau_ii(:,:,ind_k,ind_ch0(i))
          call tau_to_nu_ba3ut(mm_nu(1,1,i),ncmpl*ndim_pbmt(me_b+1),
     &                         ncmpl*ndim_pbmt(me_b+1),1)
          call tau_to_nu_ba3(mi_nu(1,1,1,i),ncmpl*nd_b_pbmt(me_b+1),
     &                       nplw_gw,ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),
     &                       1)
          call tau_to_nu_ba3ut(ii_nu(1,1,i),nn2,nn2,1)
        enddo
      endif
c ------------------------------------------------------------------
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(nn==n_pbt) then
          if(ncmpl==1) then
            call unpack_hermit3rOLD(tmp,mm_nu(1,ind_nu,1),
     &                              mi_nu(1,1,ind_nu,1),
     &                              ii_nu(1,ind_nu,1),nplwgw(k),k)
          else
            call unpack_hermit3(tmp,mm_nu(1,ind_nu,1),
     &                          mi_nu(1,1,ind_nu,1),ii_nu(1,ind_nu,1),
     &                          nplwgw(k))
          endif
          chi0=tmp
          if(ch=='___XY') then
            if(ncmpl==1) then
              call unpack_hermit3rOLD(tmp,mm_nu(1,ind_nu,2),
     &                                mi_nu(1,1,ind_nu,2),
     &                                ii_nu(1,ind_nu,2),nplwgw(k),k)
            else
              call unpack_hermit3(tmp,mm_nu(1,ind_nu,2),
     &                            mi_nu(1,1,ind_nu,2),ii_nu(1,ind_nu,2),
     &                            nplwgw(k))
            endif
            chi0=chi0+(0.d0,1.d0)*tmp
          endif
        else if(nn>n_pbt) then
          if(ncmpl==1) then
            call unpack_hermit3rOLD(tmp,mm_nu(1,ind_nu,1),
     &                              mi_nu(1,1,ind_nu,1),
     &                              ii_nu(1,ind_nu,1),nplwgw(k),k)
          else
            call unpack_hermit3(tmp,mm_nu(1,ind_nu,1),
     &                          mi_nu(1,1,ind_nu,1),ii_nu(1,ind_nu,1),
     &                          nplwgw(k))
          endif
          do j=1,n_pbt
            do i=1,n_pbt
              chi0(i,j)=tmp(i,j)
              chi0(n_pbt+i,n_pbt+j)=tmp(i,j)
            enddo
          enddo
          if(ncmpl==1) then
            call unpack_hermit3rOLD(tmp,mm_nu(1,ind_nu,2),
     &                              mi_nu(1,1,ind_nu,2),
     &                              ii_nu(1,ind_nu,2),nplwgw(k),k)
          else
            call unpack_hermit3(tmp,mm_nu(1,ind_nu,2),
     &                          mi_nu(1,1,ind_nu,2),ii_nu(1,ind_nu,2),
     &                          nplwgw(k))
          endif
          do j=1,n_pbt
            do i=1,n_pbt
              chi0(i,n_pbt+j)=tmp(i,j)
              chi0(n_pbt+i,j)=conjg(tmp(j,i))
            enddo
          enddo
        endif
c  ------- Enhanced Susceptibility ---------------------------------
        sing=.false.
        if(k==1.and.(ch=='___00'.or.ch=='00_ZZ')) then
          if(metal.and.i_nu/=0) sing=.true.
          if(.not.metal) sing=.true.
        endif
        if(sing) then
          call ixc_00_zz(ind_nu,nn,chi0,ixc0,chi)
        else
          if(k==1.and.(ch=='___00'.or.ch=='00_ZZ')) then
            ixc0(1:n_pbt,1:n_pbt)=v_tmp
          endif
          call zgemm('n','n',nn,nn,nn,(-1.d0,0.d0),chi0,nn,ixc0,nn,
     &	             (0.d0,0.d0),work,nn)
          do i=1,nn
            work(i,i)=(1.d0,0.d0)+work(i,i)
          enddo
          chi=chi0
          call zgesv(nn,nn,work,nn,ipiv,chi,nn,info)
        endif
        if(ch=='___ZZ'.or.ch=='___XY') then
          tmp1=chi(1:n_pbmt,1:n_pbmt)
          call mt_integral_pb_matrix(mtp,tmp1,k,chi_mt(1,ind_nu))
          call mt_integral_neighbors(mtp,tmp1,k,chi_neig(1,1,ind_nu))
          tmp1=chi0(1:n_pbmt,1:n_pbmt)
          call mt_integral_pb_matrix(mtp,tmp1,k,chi0_mt(1,ind_nu))
        else if(ch=='00_ZZ') then
          tmp1=chi(n_pbt+1:n_pbt+n_pbmt,n_pbt+1:n_pbt+n_pbmt)
          call mt_integral_pb_matrix(mtp,tmp1,k,chi_mt(1,ind_nu))
          call mt_integral_neighbors(mtp,tmp1,k,chi_neig(1,1,ind_nu))
          tmp1=chi0(n_pbt+1:n_pbt+n_pbmt,n_pbt+1:n_pbt+n_pbmt)
          call mt_integral_pb_matrix(mtp,tmp1,k,chi0_mt(1,ind_nu))
        endif
c ------------------------------------------------------------------
        if(iq_l/=0) then    
          call pw_from_pb(ind_k,n_pbt,nn,chi0(1,1),xi0_nu(i_nu,iq_l,1),
     &                    3)
          call pw_from_pb(ind_k,n_pbt,nn,chi(1,1),xi_nu(i_nu,iq_l,1),3)
c ----------- Inverse Macroscopic Dielectric function ----------------
          if(ch=='___00'.or.ch=='00_ZZ') then
            if(k==1) then
              xi0_nu(i_nu,iq_l,1)=0.d0
              xi_nu(i_nu,iq_l,1)=0.d0
              if(metal.and.i_nu==0) then
                eps1(i_nu,iq_l)=0.d0
                do i=1,3
                  eps_tns(i,i,i_nu)=1.d10
                enddo
              else
                p00=0.d0
                ets=pi8*p_head_rsp(:,:,ind_nu)
                do i=1,3
                  p00=p00+p_head_rsp(i,i,ind_nu)
                  ets(i,i)=(1.d0,0.d0)+ets(i,i)
                enddo
                p00=p00/3.d0
                eps1(i_nu,iq_l)=1.d0+pi8*p00
              endif
            else
              call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,
     &                   n_pbt,chi,nn,(0.d0,0.d0),work,nn)
	          do i=1,n_pbt
	            work(i,i)=(1.d0,0.d0)+work(i,i)
	          enddo
              call pw_from_pb(ind_k,n_pbt,nn,work,eps1(i_nu,iq_l),2)
            endif
          endif
          if(n_ch==3) then
            tmp=chi0(nn/2+1:nn,nn/2+1:nn)
            call pw_from_pb(ind_k,n_pbt,n_pbt,tmp,xi0_nu(i_nu,iq_l,2),3)
            tmp=chi(nn/2+1:nn,nn/2+1:nn)
            call pw_from_pb(ind_k,n_pbt,n_pbt,tmp,xi_nu(i_nu,iq_l,2),3)
            tmp=chi0(1:nn/2,nn/2+1:nn)
            call pw_from_pb(ind_k,n_pbt,n_pbt,tmp,xi0_nu(i_nu,iq_l,3),3)
            tmp=chi(1:nn/2,nn/2+1:nn)
            call pw_from_pb(ind_k,n_pbt,n_pbt,tmp,xi_nu(i_nu,iq_l,3),3)
          endif
        endif
c ------------------------------------------------------------------
        if(i_nu==0) then    
          if(ch=='___00'.or.ch=='___ZZ'.or.ch=='___XY') then
            call pw_from_pb(ind_k,n_pbt,n_pbt,chi(1,1),cc,3)
            resp_stat(k,1)=cc
          else if(ch=='00_ZZ') then
            tmp=chi(1:nn/2,1:nn/2)
            call pw_from_pb(ind_k,n_pbt,n_pbt,tmp,cc,3)
            resp_stat(k,1)=cc
            tmp=chi(nn/2+1:nn,nn/2+1:nn)
            call pw_from_pb(ind_k,n_pbt,n_pbt,tmp,cc,3)
            resp_stat(k,2)=cc
          endif
        endif
c ------------------------------------------------------------------
      enddo  !! over ind_nu
      deallocate(mm_nu,mi_nu,ii_nu,a,d,ipiv,work,ixc0,chi0,tmp,tmp1,chi)
      if(ch=='___00'.or.ch=='00_ZZ') deallocate(v_tmp)
      end
