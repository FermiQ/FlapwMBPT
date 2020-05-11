      subroutine ws_k_nu_solid_0(ind_k,w_last,w_red)
	  use atom_mod
	  use etot_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use scalapack_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer, intent(in) :: ind_k
	  logical, intent(in) :: w_last,w_red
#ifdef SCLP
      integer :: ii,jj,j
      real*8, allocatable :: vlocr(:,:)
      complex*16, allocatable :: vloc(:,:)
#endif
      logical :: pvr
      integer :: k,n_pbt,i,ind_nu,i_nu,k0,info,n,nn2,nrow,ncol
      real*8 :: ccr
	  complex*16 :: cc
	  integer, allocatable :: ipiv(:)
      real*8, allocatable :: v_tmpr(:,:),workr(:,:),v_pr(:,:),w_sr(:,:)
      complex*16, allocatable :: v_tmp(:,:),v_p(:,:),w_s(:,:),work(:,:),
     &                           pvrt(:,:),w_b(:,:)
      nn2=ncmpl*nd_b_pbint(me_b+1)
      pvr=.false.
      if(ubi=='psi'.and.psi_p/='00000') pvr=.true.
      if(ubi=='bsp'.and.bse_kernel_p/='0000') pvr=.true.
      if(ubi==' qp'.and.w_sc_qp=='non'.and.ubi_0==' gw'
     &     .and.w_sc_gw=='scf'.and.bse_kernel_p/='0000') pvr=.true.
	  k=n3_mpi_k(me_k+1)+ind_k
	  k0=k_npnt_in_line(k)
	  n_pbt=n_pbmt+nplwgw(k)
      if(ncmpl==1) then
        allocate(v_tmpr(n_pbt,n_pbt))
        call v_coul_full_r(nplwgw(k),v_tmpr,ind_k,0)
      else
        allocate(v_tmp(n_pbt,n_pbt))
        call v_coul_full(nplwgw(k),v_tmp,ind_k,0)
      endif
      if(lsclp) then
#ifdef SCLP
        ncol=ncol_pw(k)
        nrow=nrow_pw(k)
        if(complex_ro) then
          allocate(vloc(nrow,ncol))
          call matrix_to_scalapack(n_pbt,v_tmp,nrow,ncol,vloc,
     &                             refrow_pw(1,k),refcol_pw(1,k))
          deallocate(v_tmp)
          allocate(v_tmp(nrow,ncol))
          v_tmp=vloc
          deallocate(vloc)
        else
          allocate(vlocr(nrow,ncol))
          call matrix_to_scalapack_r(n_pbt,v_tmpr,nrow,ncol,vlocr,
     &                               refrow_pw(1,k),refcol_pw(1,k))
          deallocate(v_tmpr)
          allocate(v_tmpr(nrow,ncol))
          v_tmpr=vlocr
          deallocate(vlocr)
        endif
#endif
      endif
c   --- Transform P_FULL from TAU to NU representation for given Q -
      if(lsclp) then
        call tau_to_nu_ba3ut(pw_sclp(1,1,1,ind_k),
     &                       ncmpl*nrow_pw_max*ncol_pw_max,
     &                       ncmpl*nrow_pw_max*ncol_pw_max,1)
      else
        call tau_to_nu_ba3ut(p_w_q_tau_mm(1,1,ind_k),
     &                       ncmpl*ndim_pbmt(me_b+1),
     &                       ncmpl*ndim_pbmt(me_b+1),1)
        call tau_to_nu_ba3(p_w_q_tau_mi(1,1,1,ind_k),
     &                     ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                     ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),1)
        call tau_to_nu_ba3ut(p_w_q_tau_ii(1,1,ind_k),nn2,nn2,1)
      endif
c ------------------------------------------------------------------
	  do ind_nu=1,ndim3_nu
	    i_nu=me_t*ndim3_nu+ind_nu-1
        if(ncmpl==1) then
          if(lsclp) then
            allocate(workr(nrow,ncol))
            call dlacpy('a',nrow,ncol,pw_sclp(1,1,ind_nu,ind_k),
     &                  nrow_pw_max,workr,nrow)
          else
            allocate(workr(n_pbt,n_pbt))
            call unpack_hermit3r(workr,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                           p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                           p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
          endif
        else
          if(lsclp) then
            allocate(work(nrow,ncol))
            call zlacpy('a',nrow,ncol,pw_sclp(1,1,ind_nu,ind_k),
     &                  nrow_pw_max,work,nrow)
          else
            allocate(work(n_pbt,n_pbt))
            call unpack_hermit3(work,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                          p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                          p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
          endif
        endif
        if(pvr) then
          if(lsclp) then
          else
            n=n_pbmt_red+nplwgw_red(k)
            allocate(pvrt(n,n))
            call unpack_hermit(pvrt,p_vrt_nu_q(1,1,ind_nu,ind_k),n,
     &                         n_pbtot_red,n,(0.d0,0.d0),(1.d0,0.d0))
            if(ncmpl==1) then
              allocate(v_p(n_pbt,n_pbt))
              call fpb_from_rpb_full(pvrt,v_p,k,n,n_pbt,0)
              allocate(v_pr(n_pbt,n_pbt))
              call pb_r_from_c(n_pbt,n_pbt,n_pbt,v_p,v_pr,pnt(1,k))
              deallocate(v_p)
              workr=workr+v_pr
              deallocate(v_pr)
            else
              call fpb_from_rpb_full(pvrt,work,k,n,n_pbt,1)
            endif
            deallocate(pvrt)
          endif
        endif
        if(lsclp) then
#ifdef SCLP
          if(ncmpl==2) then
            if(k0/=0) call pw_from_pb_sclp(ind_k,nrow,ncol,
     &                                     refrow_pw(1,k),
     &                                     refcol_pw(1,k),work,cc,3)
          else
            if(k0/=0) call pw_from_pb_sclp_r(ind_k,nrow,ncol,
     &                                       refrow_pw(1,k),
     &                                       refcol_pw(1,k),workr,ccr,3)
          endif
c  ------- V_P = V*P -----------------------------------------------
          if(ncmpl==2) then
            allocate(v_p(nrow,ncol))
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,1,1,
     &                  desc_pw(1,k),work,1,1,desc_pw(1,k),(0.d0,0.d0),
     &                  v_p,1,1,desc_pw(1,k))
          else
            allocate(v_pr(nrow,ncol))
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,1.d0,v_tmpr,1,1,
     &                  desc_pw(1,k),workr,1,1,desc_pw(1,k),0.d0,v_pr,1,
     &                  1,desc_pw(1,k))
          endif
          if(ncmpl==2) then
            if(k0/=0) p_output(i_nu,k0)=real(cc)
c --------- V_P = 1 - V*P ------------------------------------------
            v_p=-v_p
            do j=1,ncol
              jj=refcol_pw(j,k)
              do i=1,nrow
                ii=refrow_pw(i,k)
                if(ii==jj) v_p(i,j)=(1.d0,0.d0)+v_p(i,j)
              enddo
            enddo
c ------- Contribution to eps for output ---------------------------
            if(k0/=0) then
              call pw_from_pb_sclp(ind_k,nrow,ncol,refrow_pw(1,k),
     &                             refcol_pw(1,k),v_p,cc,2)
              e_output(i_nu,k0)=real(cc)
            endif
c --------- W_s = (1 - V*P)^(-1) * V   -----------------------------
            allocate(w_s(nrow,ncol))
            w_s=v_tmp
            allocate(ipiv(nrow+nblock_pw))
            call pzgesv(n_pbt,n_pbt,v_p,1,1,desc_pw(1,k),ipiv,w_s,1,1,
     &                  desc_pw(1,k),info)
            deallocate(v_p,ipiv)
          else
            if(k0/=0) p_output(i_nu,k0)=ccr
c --------- V_P = 1 - V*P ------------------------------------------
            v_pr=-v_pr
            do j=1,ncol
              jj=refcol_pw(j,k)
              do i=1,nrow
                ii=refrow_pw(i,k)
                if(ii==jj) v_pr(i,j)=1.d0+v_pr(i,j)
              enddo
            enddo
c ------- Contribution to eps for output ---------------------------
            if(k0/=0) then
              call pw_from_pb_sclp_r(ind_k,nrow,ncol,refrow_pw(1,k),
     &                               refcol_pw(1,k),v_pr,ccr,2)
              e_output(i_nu,k0)=ccr
            endif
c --------- W_s = (1 - V*P)^(-1) * V   -----------------------------
            allocate(w_sr(nrow,ncol))
            w_sr=v_tmpr
            allocate(ipiv(nrow+nblock_pw))
            call pdgesv(n_pbt,n_pbt,v_pr,1,1,desc_pw(1,k),ipiv,w_sr,1,1,
     &                  desc_pw(1,k),info)
            deallocate(v_pr,ipiv)
          endif
          if(ncmpl==2) then
            w_s=w_s-v_tmp
c ------- Contribution to W_S for output ---------------------------
            if(k0/=0) then
              call pw_from_pb_sclp(ind_k,nrow,ncol,refrow_pw(1,k),
     &                             refcol_pw(1,k),w_s,cc,1)
              w_output(i_nu,k0)=real(cc)
            endif
            if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last)) then
c ----- Contribution to Etot_C from PW -----------------------------
              call etot_c_from_pw_sclp(nrow,ncol,k,i_nu,work,w_s)
c ----- Contribution to Phi_c in RPA  ------------------------------
              call phi_rpa_c_sclp(n_pbt,nrow,ncol,k,i_nu,work,v_tmp,
     &                            desc_pw(1,k),refrow_pw(1,k),
     &                            refcol_pw(1,k))
            endif
            work=w_s
            deallocate(w_s)
          else
            w_sr=w_sr-v_tmpr
c ------- Contribution to W_S for output ---------------------------
            if(k0/=0) then
              call pw_from_pb_sclp_r(ind_k,nrow,ncol,refrow_pw(1,k),
     &                               refcol_pw(1,k),w_sr,ccr,1)
              w_output(i_nu,k0)=ccr
            endif
            if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last)) then
c ----- Contribution to Etot_C from PW -----------------------------
              call etot_c_from_pw_sclp_r(nrow,ncol,k,i_nu,workr,w_sr)
c ----- Contribution to Phi_c in RPA  ------------------------------
              call phi_rpa_c_sclp_r(n_pbt,nrow,ncol,k,i_nu,workr,v_tmpr,
     &                              desc_pw(1,k),refrow_pw(1,k),
     &                              refcol_pw(1,k))
            endif
            workr=w_sr
            deallocate(w_sr)
          endif
#endif
        else
          if(ncmpl==2) then
            if(k0/=0) call pw_from_pb(ind_k,n_pbt,n_pbt,work,cc,3)
c  ------- V_P = V*P -----------------------------------------------
            allocate(v_p(n_pbt,n_pbt))
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,
     &                 n_pbt,work,n_pbt,(0.d0,0.d0),v_p,n_pbt)
          else
            if(k0/=0) call pw_from_pb_r(ind_k,n_pbt,n_pbt,workr,ccr,3)
c  ------- V_P = V*P -----------------------------------------------
            allocate(v_pr(n_pbt,n_pbt))
            call dgemm('n','n',n_pbt,n_pbt,n_pbt,1.d0,v_tmpr,n_pbt,
     &                 workr,n_pbt,0.d0,v_pr,n_pbt)
          endif
          if(ncmpl==2) then
            if(k0/=0) p_output(i_nu,k0)=real(cc)
c --------- V_P = 1 - V*P ------------------------------------------
            v_p=-v_p
	        do i=1,n_pbt
	          v_p(i,i)=(1.d0,0.d0)+v_p(i,i)
	        enddo
c ------- Contribution to eps for output ---------------------------
            if(k0/=0) then
              call pw_from_pb(ind_k,n_pbt,n_pbt,v_p,cc,2)
              e_output(i_nu,k0)=real(cc)
            endif
c --------- W_s = (1 - V*P)^(-1) * V   -----------------------------
            allocate(w_s(n_pbt,n_pbt))
            w_s=v_tmp
            allocate(ipiv(n_pbt))
            call zgesv(n_pbt,n_pbt,v_p,n_pbt,ipiv,w_s,n_pbt,info)
            deallocate(ipiv,v_p)
          else
            if(k0/=0) p_output(i_nu,k0)=ccr
c --------- V_P = 1 - V*P ------------------------------------------
            v_pr=-v_pr
            do i=1,n_pbt
              v_pr(i,i)=1.d0+v_pr(i,i)
            enddo
c ------- Contribution to eps for output ---------------------------
            if(k0/=0) then
              call pw_from_pb_r(ind_k,n_pbt,n_pbt,v_pr,ccr,2)
              e_output(i_nu,k0)=ccr
            endif
c --------- W_s = (1 - V*P)^(-1) * V   -----------------------------
            allocate(w_sr(n_pbt,n_pbt))
            w_sr=v_tmpr
            allocate(ipiv(n_pbt))
            call dgesv(n_pbt,n_pbt,v_pr,n_pbt,ipiv,w_sr,n_pbt,info)
            deallocate(ipiv,v_pr)
          endif
          if(ncmpl==2) then
            w_s=w_s-v_tmp
          if(k==2) then
              if(i_nu==1) then
              do i=n_pbmt+1,n_pbt
              print *,i,real(v_tmp(i,i)),real(w_s(i,i)),'k=2'
              enddo
              endif
          endif
c ------- Contribution to W_S for output ---------------------------
            if(k0/=0) then
              call pw_from_pb(ind_k,n_pbt,n_pbt,w_s,cc,1)
              w_output(i_nu,k0)=real(cc)
            endif
	        if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last)) then
c ----- Contribution to Etot_C from PW --------------------------------
	          call etot_c_from_pw(n_pbt,k,i_nu,work,w_s)
c ----- Contribution to Phi_c in RPA  ------------------------------
              call phi_rpa_c(n_pbt,n_pbt,ind_k,i_nu,work,v_tmp)
            endif
            work=w_s
            deallocate(w_s)
          else
            w_sr=w_sr-v_tmpr
          if(k==2) then
              if(i_nu==1) then
              do i=n_pbmt+1,n_pbt
              print *,i,v_tmpr(i,i),w_sr(i,i),'k=2'
              enddo
              endif
          endif
c ------- Contribution to W_S for output ---------------------------
            if(k0/=0) then
              call pw_from_pb_r(ind_k,n_pbt,n_pbt,w_sr,ccr,1)
              w_output(i_nu,k0)=ccr
            endif
            if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last)) then
c ----- Contribution to Etot_C from PW --------------------------------
              call etot_c_from_pw_r(n_pbt,k,i_nu,workr,w_sr)
c ----- Contribution to Phi_c in RPA  ------------------------------
              call phi_rpa_c_r(n_pbt,n_pbt,ind_k,i_nu,workr,v_tmpr)
            endif
            workr=w_sr
            deallocate(w_sr)
          endif
        endif
c ------------------------------------------------------------------
        if(w_red) then
          if(lsclp) then
          else
            if(ncmpl==1) then
              allocate(work(n_pbt,n_pbt))
              call pb_c_from_r(n_pbt,n_pbt,n_pbt,work,workr,pnt(1,k))
            endif
            allocate(w_b(n_pbtot_red,n_pbtot_red))
            call rpb_from_fpb_full(w_b,work,k,n_pbtot_red,n_pbt)
            if(ncmpl==1) deallocate(work)
            n=n_pbmt_red+nplwgw_red(k)
            call pack_hermit(w_b,w_nu_q_b(1,1,ind_nu,ind_k),n,
     &                       n_pbtot_red,n_pbtot_red,0.d0,1.d0)
            deallocate(w_b)
          endif
        endif
        if(lsclp) then
          if(ncmpl==1) then
            call dlacpy('a',nrow,ncol,workr,nrow,
     &                  pw_sclp(1,1,ind_nu,ind_k),nrow_pw_max)
            deallocate(workr)
          else
            call zlacpy('a',nrow,ncol,work,nrow,
     &                  pw_sclp(1,1,ind_nu,ind_k),nrow_pw_max)
            deallocate(work)
          endif
        else
          if(ncmpl==1) then
            call pack_hermit3r(workr,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                         p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                         p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
            deallocate(workr)
          else
            call pack_hermit3(work,p_w_q_tau_mm(1,ind_nu,ind_k),
     &                        p_w_q_tau_mi(1,1,ind_nu,ind_k),
     &                        p_w_q_tau_ii(1,ind_nu,ind_k),nplwgw(k))
            deallocate(work)
          endif
        endif
c ------------------------------------------------------------------
	  enddo  !! over ind_nu
      if(ncmpl==1) deallocate(v_tmpr)
	  if(ncmpl==2) deallocate(v_tmp)
c   ------- Transform NU ---> TAU representation -------------------
      if(lsclp) then
        call nu_to_tau_ba3ut(pw_sclp(1,1,1,ind_k),
     &                       ncmpl*nrow_pw_max*ncol_pw_max,
     &                       ncmpl*nrow_pw_max*ncol_pw_max)
      else
        call nu_to_tau_ba3ut(p_w_q_tau_mm(1,1,ind_k),
     &                       ncmpl*ndim_pbmt(me_b+1),
     &                       ncmpl*ndim_pbmt(me_b+1))
        call nu_to_tau_ba3(p_w_q_tau_mi(1,1,1,ind_k),
     &                     ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                     ncmpl*nd_b_pbmt(me_b+1),nplwgw(k))
        call nu_to_tau_ba3ut(p_w_q_tau_ii(1,1,ind_k),nn2,nn2)
      endif
      end
