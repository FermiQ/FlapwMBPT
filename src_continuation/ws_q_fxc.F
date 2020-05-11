      subroutine ws_q_fxc(ind_k)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_k
      integer :: k,n_pbt,i,ind_nu,i_nu,k0,info,nn2,j
      real*8 :: ccr,pi2a,v00,p0r,r0
      complex*16 :: cc
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: v_tmpr(:,:),workr(:,:),v_pr(:,:),w_sr(:,:),
     &                       fxcr(:,:),fxc0r(:,:)
      complex*16, allocatable :: v_tmp(:,:),v_p(:,:),w_s(:,:),work(:,:),
     &                           fxc(:,:),fxc0(:,:)
      pi2a=(pi+pi)/par
      nn2=ncmpl*nd_b_pbint(me_b+1)
      k=n3_mpi_k(me_k+1)+ind_k
      k0=k_npnt_in_line(k)
      n_pbt=n_pbmt+nplwgw(k)
      v00=dot_product(pnt(:,k),pnt(:,k))*pi2a**2
      v00=8.d0*pi/v00
      if(ncmpl==1) then
        allocate(v_tmpr(n_pbt,n_pbt))
        call v_coul_full_r(nplwgw(k),v_tmpr,ind_k,0)
      else
        allocate(v_tmp(n_pbt,n_pbt))
        call v_coul_full(nplwgw(k),v_tmp,ind_k,0)
      endif
c   --- Transform P_FULL from TAU to NU representation for given Q -
      call tau_to_nu_ba3ut(p_w_q_tau_mm(1,1,ind_k),
     &                     ncmpl*ndim_pbmt(me_b+1),
     &                     ncmpl*ndim_pbmt(me_b+1),1)
      call tau_to_nu_ba3(p_w_q_tau_mi(1,1,1,ind_k),
     &                   ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                   ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),1)
      call tau_to_nu_ba3ut(p_w_q_tau_ii(1,1,ind_k),nn2,nn2,1)
      if(ncmpl==1) then
        allocate(fxc0r(n_pbt,n_pbt))
        call unpack_fxc_r(n_pbt,fxc0r,fxc_bts(1,ind_k))
        fxc0r=fxc0r+v_tmpr
        allocate(fxcr(n_pbt,n_pbt))
        fxcr=0.d0
      else
        allocate(fxc0(n_pbt,n_pbt))
        call unpack_fxc(n_pbt,fxc0,fxc_bts(1,ind_k))
        fxc0=fxc0+v_tmp
        allocate(fxc(n_pbt,n_pbt))
        fxc=(0.d0,0.d0)
      endif
c ------------------------------------------------------------------
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(ncmpl==1) then
          allocate(workr(n_pbt,n_pbt))
          call unpack_hermit3r(workr,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                         p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                         p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
        else
          allocate(work(n_pbt,n_pbt))
          call unpack_hermit3(work,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                        p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                        p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
        endif
        if(ncmpl==2) then
          call pw_from_pb(ind_k,n_pbt,n_pbt,work,cc,3)
          p0r=cc
c  ------- V_P = V*P -----------------------------------------------
          allocate(v_p(n_pbt,n_pbt))
          call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),work,
     &               n_pbt,fxc0,n_pbt,(0.d0,0.d0),v_p,n_pbt)
        else
          call pw_from_pb_r(ind_k,n_pbt,n_pbt,workr,ccr,3)
          p0r=ccr
c  ------- V_P = V*P -----------------------------------------------
          allocate(v_pr(n_pbt,n_pbt))
          call dgemm('n','n',n_pbt,n_pbt,n_pbt,1.d0,workr,n_pbt,
     &               fxc0r,n_pbt,0.d0,v_pr,n_pbt)
        endif
        if(ncmpl==2) then
          if(k0/=0) p_output(i_nu,k0)=real(cc)
c --------- V_P = 1 - V*P ------------------------------------------
          v_p=-v_p
          do i=1,n_pbt
            v_p(i,i)=(1.d0,0.d0)+v_p(i,i)
          enddo
c --------- R = (1 - R0*(V+f))^(-1) * R0   -------------------------
          allocate(ipiv(n_pbt))
          call zgesv(n_pbt,n_pbt,v_p,n_pbt,ipiv,work,n_pbt,info)
          deallocate(ipiv,v_p)
c --------- W_s = V * R * V   --------------------------------------
          allocate(w_s(n_pbt,n_pbt))
          call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,n_pbt,
     &               work,n_pbt,(0.d0,0.d0),w_s,n_pbt)
c ------- Contribution to eps^(-1) for output ----------------------
          if(k0/=0) then
            call pw_from_pb(ind_k,n_pbt,n_pbt,w_s,cc,2)
            e_output(i_nu,k0)=1.d0+real(cc)
          endif
          work=w_s
          call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),work,n_pbt,
     &               v_tmp,n_pbt,(0.d0,0.d0),w_s,n_pbt)
          deallocate(work)
c ------- Contribution to W_S for output ---------------------------
          if(k0/=0) then
            call pw_from_pb(ind_k,n_pbt,n_pbt,w_s,cc,1)
            w_output(i_nu,k0)=real(cc)
          endif
        else
          if(k0/=0) p_output(i_nu,k0)=ccr
c --------- V_P = 1 - V*P ------------------------------------------
          v_pr=-v_pr
          do i=1,n_pbt
            v_pr(i,i)=1.d0+v_pr(i,i)
          enddo
c --------- R = (1 - R0*(V+f))^(-1) * R0   -------------------------
          allocate(ipiv(n_pbt))
          call dgesv(n_pbt,n_pbt,v_pr,n_pbt,ipiv,workr,n_pbt,info)
          deallocate(ipiv,v_pr)
c --------- W_s = V * R * V   --------------------------------------
          allocate(w_sr(n_pbt,n_pbt))
          call dgemm('n','n',n_pbt,n_pbt,n_pbt,1.d0,v_tmpr,n_pbt,workr,
     &               n_pbt,0.d0,w_sr,n_pbt)
c ------- Contribution to eps^(-1) for output ----------------------
          if(k0/=0) then
            call pw_from_pb_r(ind_k,n_pbt,n_pbt,w_sr,ccr,2)
            e_output(i_nu,k0)=1.d0+ccr
          endif
          workr=w_sr
          call dgemm('n','n',n_pbt,n_pbt,n_pbt,1.d0,workr,n_pbt,v_tmpr,
     &               n_pbt,0.d0,w_sr,n_pbt)
          deallocate(workr)
c ------- Contribution to W_S for output ---------------------------
          if(k0/=0) then
            call pw_from_pb_r(ind_k,n_pbt,n_pbt,w_sr,ccr,1)
            w_output(i_nu,k0)=ccr
          endif
        endif
c ------------------------------------------------------------------
        if(ncmpl==1) then
          call pack_hermit3r(w_sr,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                       p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                       p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
          if(i_nu==0) fxcr=(v_tmpr+w_sr)/v00/p0r
          deallocate(w_sr)
        else
          call pack_hermit3(w_s,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                      p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                      p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
          if(i_nu==0) fxc=(v_tmp+w_s)/v00/p0r
          deallocate(w_s)
        endif
c ------------------------------------------------------------------
      enddo  !! over ind_nu
      if(nproc_t*nproc_b/=1) then
        if(ncmpl==1) then
          call brdcst(msgdbl,fxcr,8*n_pbt**2,0,comm_t_b)
        else
          call brdcst(msgdbl,fxc,16*n_pbt**2,0,comm_t_b)
        endif
      endif
c ------------------ Mixing ----------------------------------------
      if(ncmpl==1) then
        do j=1,n_pbt
          do i=1,n_pbt
            r0=fxc0r(i,j)-v_tmpr(i,j)
            fxcr(i,j)=0.3d0*fxcr(i,j)+0.7d0*r0
          enddo
        enddo
      else
        do j=1,n_pbt
          do i=1,n_pbt
            cc=fxc0(i,j)-v_tmp(i,j)
            fxc(i,j)=0.3d0*fxc(i,j)+0.7d0*cc
          enddo
        enddo
      endif
c ------------------------------------------------------------------
      if(ncmpl==1) then
        call pack_fxc_r(n_pbt,fxcr,fxc_bts(1,ind_k))
      else
        call pack_fxc(n_pbt,fxc,fxc_bts(1,ind_k))
      endif
      if(ncmpl==1) deallocate(v_tmpr,fxcr,fxc0r)
      if(ncmpl==2) deallocate(v_tmp,fxc,fxc0)
c   ------- Transform NU ---> TAU representation -------------------
      call nu_to_tau_ba3ut(p_w_q_tau_mm(1,1,ind_k),
     &                     ncmpl*ndim_pbmt(me_b+1),
     &                     ncmpl*ndim_pbmt(me_b+1))
      call nu_to_tau_ba3(p_w_q_tau_mi(1,1,1,ind_k),
     &                   ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                   ncmpl*nd_b_pbmt(me_b+1),nplwgw(k))
      call nu_to_tau_ba3ut(p_w_q_tau_ii(1,1,ind_k),nn2,nn2)
      end
