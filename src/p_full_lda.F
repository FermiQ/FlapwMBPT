      subroutine p_full_lda(ch,ind_k,n_ch0,ind_ch0,n_ix0,ind_ix0)
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
     &                       ind_ix0(n_ix0)
      integer :: k,n_pbt,i,ind_nu,i_nu,info,n,nn2,i00,i0z,izz
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: tmpr(:,:),thet_q(:,:,:),r00r(:,:),
     &                       r0zr(:,:),a11r(:,:),a12r(:,:),a21r(:,:),
     &                       a22r(:,:)
      complex*16, allocatable :: tmp(:,:),r00(:,:),r0z(:,:),a11(:,:),
     &                           a12(:,:),a21(:,:),a22(:,:)
      k=n3_mpi_k(me_k+1)+ind_k
      n_pbt=n_pbmt+nplwgw(k)
      nn2=ncmpl*nd_b_pbint(me_b+1)
      if(ch=='___00') n=1
      if(ch=='00_ZZ') n=2
      allocate(ipiv(n_pbt))
      i00=ind_ixc(1)
      if(n==2) then
        i0z=ind_ixc(3)
        izz=ind_ixc(2)
      endif
c ------------------------------------------------------------------
      allocate(thet_q(ncmpl*n_pbt,n_pbt,n_ixc_0))
      call theta_q_ss_0(k,n_pbt,thet_q)
      if(k==1) call timel('**** F-xc finished *****************')
c   --- Transform P_FULL from TAU to NU representation for given Q -
      call tau_to_nu_ba3ut(p_w_q_tau_mm(1,1,ind_k),
     &                     ncmpl*ndim_pbmt(me_b+1),
     &                     ncmpl*ndim_pbmt(me_b+1),1)
      call tau_to_nu_ba3(p_w_q_tau_mi(1,1,1,ind_k),
     &                   ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                   ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),1)
      call tau_to_nu_ba3ut(p_w_q_tau_ii(1,1,ind_k),nn2,nn2,1)
      if(n==2) then
        call tau_to_nu_ba3ut(p_0z_mm(1,1,ind_k),ncmpl*ndim_pbmt(me_b+1),
     &                       ncmpl*ndim_pbmt(me_b+1),1)
        call tau_to_nu_ba3(p_0z_mi(1,1,1,ind_k),
     &                     ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                     ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),1)
        call tau_to_nu_ba3ut(p_0z_ii(1,1,ind_k),nn2,nn2,1)
      endif
c ------------------------------------------------------------------
      if(k==1) call timel('**** Nu-loop starts ****************')
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(ncmpl==1) then
          allocate(r00r(n_pbt,n_pbt))
          call unpack_hermit3r(r00r,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                         p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                         p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
          allocate(a11r(n_pbt,n_pbt))
          call dgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r00r,n_pbt,
     &               thet_q(1,1,i00),n_pbt,0.d0,a11r,n_pbt)
          do i=1,n_pbt
            a11r(i,i)=1.d0+a11r(i,i)
          enddo
          if(n==1) then
            call dgesv(n_pbt,n_pbt,a11r,n_pbt,ipiv,r00r,n_pbt,info)
            deallocate(a11r)
          else
            allocate(r0zr(n_pbt,n_pbt))
            call unpack_hermit3r(r0zr,p_0z_mm(1,ind_nu,ind_k),
     &                           p_0z_mi(1,1,ind_nu,ind_k),
     &                           p_0z_ii(1,ind_nu,ind_k),nplwgw(k))
            allocate(a22r(n_pbt,n_pbt))
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r0zr,n_pbt,
     &                 thet_q(1,1,i0z),n_pbt,0.d0,a22r,n_pbt)
            a11r=a11r+a22r
            do i=1,n_pbt
              a22r(i,i)=1.d0+a22r(i,i)
            enddo
            allocate(a12r(n_pbt,n_pbt))
            allocate(a21r(n_pbt,n_pbt))
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r00r,n_pbt,
     &                 thet_q(1,1,i0z),n_pbt,0.d0,a12r,n_pbt)
            a21r=a12r
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r0zr,n_pbt,
     &                 thet_q(1,1,izz),n_pbt,1.d0,a12r,n_pbt)
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r0zr,n_pbt,
     &                 thet_q(1,1,i00),n_pbt,1.d0,a21r,n_pbt)
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r00r,n_pbt,
     &                 thet_q(1,1,izz),n_pbt,1.d0,a22r,n_pbt)
            call invers_r(n_pbt,a22r,n_pbt)
            allocate(tmpr(n_pbt,n_pbt))
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,1.d0,a12r,n_pbt,a22r,
     &                 n_pbt,0.d0,tmpr,n_pbt)
            deallocate(a12r,a22r)
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,tmpr,n_pbt,a21r,
     &                 n_pbt,1.d0,a11r,n_pbt)
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,tmpr,n_pbt,r0zr,
     &                 n_pbt,1.d0,r00r,n_pbt)
            deallocate(tmpr,a21r,r0zr)
            call dgesv(n_pbt,n_pbt,a11r,n_pbt,ipiv,r00r,n_pbt,info)
            deallocate(a11r)
          endif
          call pack_hermit3r(r00r,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                       p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                       p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
          deallocate(r00r)
        else if(ncmpl==2) then
          allocate(r00(n_pbt,n_pbt))
          call unpack_hermit3(r00,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                        p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                        p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
          allocate(a11(n_pbt,n_pbt))
          call zgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r00,n_pbt,
     &               thet_q(1,1,i00),n_pbt,(0.d0,0.d0),a11,n_pbt)
          do i=1,n_pbt
            a11(i,i)=(1.d0,0.d0)+a11(i,i)
          enddo
          if(n==1) then
            call zgesv(n_pbt,n_pbt,a11,n_pbt,ipiv,r00,n_pbt,info)
            deallocate(a11)
          else
            allocate(r0z(n_pbt,n_pbt))
            call unpack_hermit3(r0z,p_0z_mm(1,ind_nu,ind_k),
     &                          p_0z_mi(1,1,ind_nu,ind_k),
     &                          p_0z_ii(1,ind_nu,ind_k),nplwgw(k))
            allocate(a22(n_pbt,n_pbt))
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r0z,n_pbt,
     &                 thet_q(1,1,i0z),n_pbt,(0.d0,0.d0),a22,n_pbt)
            a11=a11+a22
            do i=1,n_pbt
              a22(i,i)=(1.d0,0.d0)+a22(i,i)
            enddo
            allocate(a12(n_pbt,n_pbt))
            allocate(a21(n_pbt,n_pbt))
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r00,n_pbt,
     &                 thet_q(1,1,i0z),n_pbt,(0.d0,0.d0),a12,n_pbt)
            a21=a12
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r0z,n_pbt,
     &                 thet_q(1,1,izz),n_pbt,(1.d0,0.d0),a12,n_pbt)
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r0z,n_pbt,
     &                 thet_q(1,1,i00),n_pbt,(1.d0,0.d0),a21,n_pbt)
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r00,n_pbt,
     &                 thet_q(1,1,izz),n_pbt,(1.d0,0.d0),a22,n_pbt)
            call invers_h(n_pbt,a22,n_pbt)
            allocate(tmp(n_pbt,n_pbt))
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),a12,n_pbt,
     &                 a22,n_pbt,(0.d0,0.d0),tmp,n_pbt)
            deallocate(a12,a22)
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),tmp,n_pbt,
     &                 a21,n_pbt,(1.d0,0.d0),a11,n_pbt)
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),tmp,n_pbt,
     &                 r0z,n_pbt,(1.d0,0.d0),r00,n_pbt)
            deallocate(tmp,a21,r0z)
            call zgesv(n_pbt,n_pbt,a11,n_pbt,ipiv,r00,n_pbt,info)
            deallocate(a11)
          endif
          call pack_hermit3(r00,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                      p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                      p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
          deallocate(r00)
        endif
      enddo  !! over ind_nu
      if(k==1) call timel('**** Nu-loop finished **************')
c   ------- Transform NU ---> TAU representation -------------------
      call nu_to_tau_ba3ut(p_w_q_tau_mm(1,1,ind_k),
     &                     ncmpl*ndim_pbmt(me_b+1),
     &                     ncmpl*ndim_pbmt(me_b+1))
      call nu_to_tau_ba3(p_w_q_tau_mi(1,1,1,ind_k),
     &                   ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                   ncmpl*nd_b_pbmt(me_b+1),nplwgw(k))
      call nu_to_tau_ba3ut(p_w_q_tau_ii(1,1,ind_k),nn2,nn2)
      deallocate(ipiv)
      if(k==1) call timel('**** P_FULL_LDA finished ***********')
      end
