      subroutine ws0_q0(w_last,w_red,eout,e1out)
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
      logical, intent(in) :: w_last,w_red
      complex*16, intent(out) :: eout(3,3,0:n_nu),e1out(3,3,0:n_nu)
#ifdef SCLP
      integer :: nrow,ncol,myrow,mycol,desc(9),j1,i1,k,numroc,info,
     &           desc_bm(9),desc_t(9)
      integer, allocatable :: refrow(:),refcol(:)
      real*8, allocatable :: bm_sclpr(:,:)
      complex*16, allocatable :: bm_sclp(:,:)
#endif
      logical :: pvr
      integer :: n_pbt,i,ind_nu,i_nu,nne,ii,j,n,nn2,idum,jj
      real*8 :: pi2a,pi8,p0,sqpi4,p00,pws,cc1r
      complex*16 :: cc,cc1
      real*8, allocatable :: workr(:,:),w_sr(:,:),tmpr(:,:),p01r(:),
     &                       p_headr(:,:,:),p_wingr(:,:,:),ww(:,:),
     &                       vsq(:)
      complex*16, allocatable :: w_s(:,:),tmp(:,:),work(:,:),
     &                           w_b(:,:),p_head(:,:,:),p_wing(:,:,:),
     &                           p01(:),pvrt(:,:)
      pi2a=(pi+pi)/par
      pi8=8.d0*pi
      sqpi4=sqrt(4.d0*pi)
      nn2=ncmpl*nd_b_pbint(me_b+1)
      pvr=.false.
      if(ubi=='psi'.and.psi_p/='00000') pvr=.true.
      if(ubi=='bsp'.and.bse_kernel_p/='0000') pvr=.true.
      if(ubi==' qp'.and.w_sc_qp=='non'.and.ubi_0==' gw'
     &     .and.w_sc_gw=='scf'.and.bse_kernel_p/='0000') pvr=.true.
      n_pbt=n_pbmt+nplwgw(1)
      nne=n_opt_pb
      if(lsclp) then
#ifdef SCLP
        call blacs_gridinfo(sctxt_pw,nrow_b,ncol_b,myrow,mycol)
        nrow=numroc(nne,nblock_pw,myrow,0,nrow_b)
        ncol=numroc(nne,nblock_pw,mycol,0,ncol_b)
        allocate(refrow(nrow))
        allocate(refcol(ncol))
        do i=1,ncol
          k=(i-1)/nblock_pw
          k=(k*ncol_b+mycol)*nblock_pw+i-k*nblock_pw
          refcol(i)=k
        enddo
        do i=1,nrow
          k=(i-1)/nblock_pw
          k=(k*nrow_b+myrow)*nblock_pw+i-k*nblock_pw
          refrow(i)=k
        enddo
        call descinit(desc,nne,nne,nblock_pw,nblock_pw,0,0,sctxt_pw,
     &                nrow,info)
        call descinit(desc_bm,nne,n_pbt,nblock_pw,nblock_pw,0,0,
     &                sctxt_pw,nrow,info)
        call descinit(desc_t,n_pbt,nne,nblock_pw,nblock_pw,0,0,
     &                sctxt_pw,nrow_pw(1),info)
        if(ncmpl==1) then
          allocate(bm_sclpr(nrow,ncol_pw(1)))
          do j=1,ncol_pw(1)
            jj=refcol_pw(j,1)
            do i=1,nrow
              ii=refrow(i)
              bm_sclpr(i,j)=bm_coef_r(ii,jj)
            enddo
          enddo
          allocate(tmpr(nrow_pw(1),ncol))
        else
          allocate(bm_sclp(nrow,ncol_pw(1)))
          do j=1,ncol_pw(1)
            jj=refcol_pw(j,1)
            do i=1,nrow
              ii=refrow(i)
              bm_sclp(i,j)=bm_coef(ii,jj)
            enddo
          enddo
          allocate(tmp(nrow_pw(1),ncol))
        endif
#endif
      endif
      allocate(vsq(nne))
      do i=1,nne
        vsq(i)=sqrt(v_opt_e(i))
      enddo
c   ------- Transform P_FULL from TAU to NU representation for Q=0 -
      if(lsclp) then
        call tau_to_nu_ba3ut(pw_sclp,ncmpl*nrow_pw_max*ncol_pw_max,
     &                       ncmpl*nrow_pw_max*ncol_pw_max,1)
      else
        call tau_to_nu_ba3ut(p_w_q_tau_mm,ncmpl*ndim_pbmt(me_b+1),
     &                       ncmpl*ndim_pbmt(me_b+1),1)
        call tau_to_nu_ba3(p_w_q_tau_mi,ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                     ncmpl*nd_b_pbmt(me_b+1),nplwgw(1),1)
        call tau_to_nu_ba3ut(p_w_q_tau_ii,nn2,nn2,1)
      endif
c ------------------------------------------------------------------
      if(lsclp) then
#ifdef SCLP
        if(ncmpl==1) then
          allocate(w_sr(nrow,ncol))
        else
          allocate(w_s(nrow,ncol))
        endif
#endif
      else
        if(ncmpl==1) then
          allocate(w_sr(nne,nne))
          allocate(tmpr(n_pbt,nne))
        else
          allocate(w_s(nne,nne))
          allocate(tmp(n_pbt,nne))
        endif
      endif
c ------------------------------------------------------------------
      ws_head_nu=0.d0
      if(ncmpl==1) then
        allocate(p_headr(3,3,ndim3_nu))
        allocate(p_wingr(3,nne,ndim3_nu))
        call dcopy(3*nne*ndim3_nu,p_wing_all,1,p_wingr,1)
        do j=1,3
          do i=1,j
            ii=ind_tens33(i,j)
            p_headr(i,j,:)=p_head_all(ii,:)
            p_headr(j,i,:)=p_head_all(ii,:)
          enddo
        enddo
        if(pvr) then
          do j=1,3
            do i=1,j
              ii=ind_tens33(i,j)
              p_headr(i,j,:)=p_headr(i,j,:)+p_head_vrt(ii,:)
              p_headr(j,i,:)=p_headr(j,i,:)+p_head_vrt(ii,:)
            enddo
          enddo
        endif
        call output_p_tensor_r(p_headr)
        if(metal) allocate(p01r(nne))
      else
        allocate(p_head(3,3,ndim3_nu))
        allocate(p_wing(3,nne,ndim3_nu))
        call zcopy(3*nne*ndim3_nu,p_wing_all,1,p_wing,1)
        do j=1,3
          do i=1,j
            ii=ind_tens33(i,j)
            p_head(i,j,:)=p_head_all(ii,:)
            p_head(j,i,:)=p_head_all(ii,:)
          enddo
        enddo
        if(pvr) then
          do j=1,3
            do i=1,j
              ii=ind_tens33(i,j)
              p_head(i,j,:)=p_head(i,j,:)+p_head_vrt(ii,:)
              p_head(j,i,:)=p_head(j,i,:)+p_head_vrt(ii,:)
            enddo
          enddo
        endif
        call output_p_tensor(p_head)
        if(metal) allocate(p01(nne))
      endif
c ------------------------------------------------------------------
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(ncmpl==1) then
          if(lsclp) then
            allocate(workr(nrow_pw(1),ncol_pw(1)))
            call dlacpy('a',nrow_pw(1),ncol_pw(1),pw_sclp(1,1,ind_nu,1),
     &                  nrow_pw_max,workr,nrow_pw(1))
          else
            allocate(workr(n_pbt,n_pbt))
            call unpack_hermit3r(workr,p_w_q_tau_mm(1,ind_nu,1),
     &                           p_w_q_tau_mi(1,1,ind_nu,1),
     &                           p_w_q_tau_ii(1,ind_nu,1),nplwgw(1))
          endif
        else
          if(lsclp) then
            allocate(work(nrow_pw(1),ncol_pw(1)))
            call zlacpy('a',nrow_pw(1),ncol_pw(1),pw_sclp(1,1,ind_nu,1),
     &                  nrow_pw_max,work,nrow_pw(1))
          else
            allocate(work(n_pbt,n_pbt))
            call unpack_hermit3(work,p_w_q_tau_mm(1,ind_nu,1),
     &                          p_w_q_tau_mi(1,1,ind_nu,1),
     &                          p_w_q_tau_ii(1,ind_nu,1),nplwgw(1))
          endif
        endif
        if(pvr) then
          if(lsclp) then
          else
            n=n_pbmt_red+nplwgw_red(1)
            allocate(pvrt(n,n))
            call unpack_hermit(pvrt,p_vrt_nu_q(1,1,ind_nu,1),n,
     &                         n_pbtot_red,n,(0.d0,0.d0),(1.d0,0.d0))
            if(ncmpl==1) then
              allocate(work(n_pbt,n_pbt))
              call fpb_from_rpb_full(pvrt,work,1,n,n_pbt,0)
              allocate(ww(n_pbt,n_pbt))
              call pb_r_from_c(n_pbt,n_pbt,n_pbt,work,ww,pnt(1,1))
              deallocate(work)
              workr=workr+ww
              deallocate(ww)
            else
              call fpb_from_rpb_full(pvrt,work,1,n,n_pbt,1)
            endif
            deallocate(pvrt)
          endif
        endif
        if(lsclp) then
          if(ncmpl==1) then
            call pw_from_pb_sclp_r(1,nrow_pw(1),ncol_pw(1),refrow_pw,
     &                             refcol_pw,workr,p0,3)
          else
            call pw_from_pb_sclp(1,nrow_pw(1),ncol_pw(1),refrow_pw,
     &                           refcol_pw,work,cc,3)
            p0=cc
          endif
        else
          if(ncmpl==1) then
            call pw_from_pb_r(1,n_pbt,n_pbt,workr,p0,3)
          else
            call pw_from_pb(1,n_pbt,n_pbt,work,cc,3)
            p0=cc
          endif
        endif
c -------- Transform P to the orthonormal basis set ----------------
        if(metal) then
          if(ncmpl==1) then
            p01r=0.d0
            if(lsclp) then
              do j=1,ncol_pw(1)
                jj=refcol_pw(j,1)
                cc1r=0.d0
                do i=1,nrow_pw(1)
                  ii=refrow_pw(i,1)
                  cc1r=cc1r+pw_pb_tild_r(ii,1)*workr(i,j)
                enddo
                do i=1,nne
                  p01r(i)=p01r(i)+cc1r*bm_coef_r(i,jj)
                enddo
              enddo
              call dgop(p01r,nne,'  +',comm_b)
            else
              do j=1,n_pbt
                cc1r=0.d0
                do i=1,n_pbt
                  cc1r=cc1r+pw_pb_tild_r(i,1)*workr(i,j)
                enddo
                do i=1,nne
                  p01r(i)=p01r(i)+cc1r*bm_coef_r(i,j)
                enddo
              enddo
            endif
          else
            p01=(0.d0,0.d0)
            if(lsclp) then
              do j=1,ncol_pw(1)
                jj=refcol_pw(j,1)
                cc1=0.d0
                do i=1,nrow_pw(1)
                  ii=refrow_pw(i,1)
                  cc1=cc1+conjg(pw_pb_tild(ii,1))*work(i,j)
                enddo
                do i=1,nne
                  p01(i)=p01(i)+cc1*conjg(bm_coef(i,jj))
                enddo
              enddo
              call dgop(p01,2*nne,'  +',comm_b)
            else
              do j=1,n_pbt
                cc1=(0.d0,0.d0)
                do i=1,n_pbt
                  cc1=cc1+conjg(pw_pb_tild(i,1))*work(i,j)
                enddo
                do i=1,nne
                  p01(i)=p01(i)+cc1*conjg(bm_coef(i,j))
                enddo
              enddo
            endif
          endif
        endif
        if(lsclp) then
#ifdef SCLP
          if(ncmpl==1) then
            call pdgemm('n','t',n_pbt,nne,n_pbt,1.d0,workr,1,1,
     &                  desc_pw(1,1),bm_sclpr,1,1,desc_bm,0.d0,tmpr,1,
     &                  1,desc_t)
            call pdgemm('n','n',nne,nne,n_pbt,1.d0,bm_sclpr,1,1,
     &                  desc_bm,tmpr,1,1,desc_t,0.d0,w_sr,1,1,desc)
            deallocate(workr)
            allocate(workr(nrow,ncol))
            if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last))
     &          workr=w_sr
          else
            call pzgemm('n','c',n_pbt,nne,n_pbt,(1.d0,0.d0),work,1,1,
     &                  desc_pw(1,1),bm_sclp,1,1,desc_bm,(0.d0,0.d0),
     &                  tmp,1,1,desc_t)
            call pzgemm('n','n',nne,nne,n_pbt,(1.d0,0.d0),bm_sclp,1,1,
     &                  desc_bm,tmp,1,1,desc_t,(0.d0,0.d0),w_s,1,1,desc)
            deallocate(work)
            allocate(work(nrow,ncol))
            if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last)) work=w_s
          endif
#endif
        else
          if(ncmpl==1) then
            call dgemm('n','t',n_pbt,nne,n_pbt,1.d0,workr,n_pbt,
     &                 bm_coef_r,nne,0.d0,tmpr,n_pbt)
            call dgemm('n','n',nne,nne,n_pbt,1.d0,bm_coef_r,nne,
     &                 tmpr,n_pbt,0.d0,w_sr,nne)
            if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last)) then
              workr(1:nne,1:nne)=w_sr
            endif
          else
            call zgemm('n','c',n_pbt,nne,n_pbt,(1.d0,0.d0),work,n_pbt,
     &                 bm_coef,nne,(0.d0,0.d0),tmp,n_pbt)
            call zgemm('n','n',nne,nne,n_pbt,(1.d0,0.d0),bm_coef,nne,
     &                 tmp,n_pbt,(0.d0,0.d0),w_s,nne)
            if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last)) then
              work(1:nne,1:nne)=w_s
            endif
          endif
        endif
        p_output(i_nu,1)=p0
c ---------- Evaluate [1-V_11^(1/2)*P_11*V_11^(1/2)]^(-1) ----------
        if(lsclp) then
#ifdef SCLP
          if(ncmpl==1) then
            do j=1,ncol
              jj=refcol(j)
              do i=1,nrow
                ii=refrow(i)
                w_sr(i,j)=-vsq(ii)*w_sr(i,j)*vsq(jj)
                if(ii==jj) w_sr(i,j)=w_sr(i,j)+1.d0
              enddo
            enddo
            call invers_r_sclp(nne,nrow,ncol,nblock_pw,desc,w_sr)
          else
            do j=1,ncol
              jj=refcol(j)
              do i=1,nrow
                ii=refrow(i)
                w_s(i,j)=-vsq(ii)*w_s(i,j)*vsq(jj)
                if(ii==jj) w_s(i,j)=w_s(i,j)+(1.d0,0.d0)
              enddo
            enddo
            call invers_z_sclp(nne,nrow,ncol,nblock_pw,desc,w_s)
          endif
#endif
        else
          if(ncmpl==1) then
            do j=1,nne
              do i=1,nne
                w_sr(i,j)=-vsq(i)*w_sr(i,j)*vsq(j)
              enddo
            enddo
            do i=1,nne
              w_sr(i,i)=1.d0+w_sr(i,i)
            enddo
            call invers_r(nne,w_sr,nne)
          else
            do j=1,nne
              do i=1,nne
                w_s(i,j)=-vsq(i)*w_s(i,j)*vsq(j)
              enddo
            enddo
            do i=1,nne
              w_s(i,i)=(1.d0,0.d0)+w_s(i,i)
            enddo
            call invers_h(nne,w_s,nne)
          endif
        endif
c ------------------------------------------------------------------
        if(metal.and.i_nu==0) then
          if(lsclp) then
#ifdef SCLP
            if(ncmpl==1) then
              call wsq0_met_r(ncol,nrow,refcol,refrow,nblock_pw,desc,
     &                        nne,lsclp,p01r,vsq,p0,w_sr,nrow,ncol)
            else
              call wsq0_met(ncol,nrow,refcol,refrow,nblock_pw,desc,nne,
     &                      lsclp,p01,vsq,p0,w_s,nrow,ncol)
            endif
#endif
          else
            if(ncmpl==1) then
              call wsq0_met_r(idum,idum,idum,idum,idum,idum,nne,lsclp,
     &                        p01r,vsq,p0,w_sr,nne,nne)
            else
              call wsq0_met(idum,idum,idum,idum,idum,idum,nne,lsclp,p01,
     &                      vsq,p0,w_s,nne,nne)
            endif
          endif
          ws_head_nu(0)=-pi8
          e_output(i_nu,1)=1.d30
          e1out(:,:,0)=(0.d0,0.d0)
          eout(1,1,0)=(1.d30,0.d0)
          eout(2,2,0)=(1.d30,0.d0)
          eout(3,3,0)=(1.d30,0.d0)
          p00=p0
c ------------------------------------------------------------------
        else
          if(lsclp) then
#ifdef SCLP
            if(ncmpl==1) then
              call wsq0_ins_r(p00,p_headr,p_wingr,lsclp,ncol,nrow,
     &                        refcol,refrow,nblock_pw,desc,nne,w_sr,
     &                        nrow,ncol,eout,e1out,ind_nu,i_nu,vsq)
            else
              call wsq0_ins(p00,p_head,p_wing,lsclp,ncol,nrow,refcol,
     &                      refrow,nblock_pw,desc,nne,w_s,nrow,ncol,
     &                      eout,e1out,ind_nu,i_nu,vsq)
            endif
#endif
          else
            if(ncmpl==1) then
              call wsq0_ins_r(p00,p_headr,p_wingr,lsclp,idum,idum,idum,
     &                        idum,idum,idum,nne,w_sr,nne,nne,eout,
     &                        e1out,ind_nu,i_nu,vsq)
            else
              call wsq0_ins(p00,p_head,p_wing,lsclp,idum,idum,idum,idum,
     &                      idum,idum,nne,w_s,nne,nne,eout,e1out,
     &                      ind_nu,i_nu,vsq)
            endif
          endif
        endif
c ------------------------------------------------------------------
        if(psi_fncl_use==2.or.(psi_fncl_use==1.and.w_last)) then
c ----- Contribution to Etot_C from PW -----------------------------
          if(lsclp) then
#ifdef SCLP
            if(ncmpl==1) then
              if(metal.and.i_nu==0) then
                call p_ws_conv_q0_met_sclp_r(nne,nrow,ncol,p01r,
     &                                       ws01_met_r,p0,
     &                                       ws_head_nu(0),workr,w_sr,
     &                                       pws,q2aver)
              else
                call p_ws_conv_q0_sclp_r(nne,nrow,ncol,
     &                                   p_wingr(1,1,ind_nu),
     &                                   ws_wing_r(1,1,ind_nu),
     &                                   p_headr(1,1,ind_nu),
     &                                   ws_head_tens_r(1,1,ind_nu),
     &                                   workr,w_sr,pws)
              endif
            else
              if(metal.and.i_nu==0) then
                call p_ws_conv_q0_met_sclp(nne,nrow,ncol,p01,
     &                                     ws01_met,p0,ws_head_nu(0),
     &                                     work,w_s,pws,q2aver)
              else
                call p_ws_conv_q0_sclp(nne,nrow,ncol,p_wing(1,1,ind_nu),
     &                                 ws_wing(1,1,ind_nu),
     &                                 p_head(1,1,ind_nu),
     &                                 ws_head_tens(1,1,ind_nu),work,
     &                                 w_s,pws)
              endif
            endif
#endif
          else
            if(ncmpl==1) then
              if(metal.and.i_nu==0) then
                call p_ws_conv_q0_met_r(nne,n_pbt,nne,p01r,ws01_met_r,
     &                                  p0,ws_head_nu(0),workr,w_sr,pws,
     &                                  q2aver)
              else
                call p_ws_conv_q0_r(nne,n_pbt,nne,p_wingr(1,1,ind_nu),
     &                              ws_wing_r(1,1,ind_nu),
     &                              p_headr(1,1,ind_nu),
     &                              ws_head_tens_r(1,1,ind_nu),workr,
     &                              w_sr,pws)
              endif
            else
              if(metal.and.i_nu==0) then
                call p_ws_conv_q0_met(nne,n_pbt,nne,p01,ws01_met,p0,
     &                                ws_head_nu(0),work,w_s,pws,q2aver)
              else
                call p_ws_conv_q0(nne,n_pbt,nne,p_wing(1,1,ind_nu),
     &                            ws_wing(1,1,ind_nu),
     &                            p_head(1,1,ind_nu),
     &                            ws_head_tens(1,1,ind_nu),work,w_s,pws)
              endif
            endif
          endif
          e_c_solid_pw=e_c_solid_pw-et_sum_nu(i_nu)*wgt(1)*pws
c ----- Contribution to Phi_c in RPA  ------------------------------
          if(lsclp) then
#ifdef SCLP
            if(ncmpl==1) then
              if(metal.and.i_nu==0) then
                call phi_rpa_c_q0_met_sclp_r(nne,workr,p00,p01r,nrow,
     &                                       ncol,refrow,refcol)
              else
                call phi_rpa_c_q0_sclp_r(nne,i_nu,workr,
     &                                   p_headr(1,1,ind_nu),
     &                                   p_wingr(1,1,ind_nu),nrow,ncol,
     &                                   refrow,refcol)
              endif
            else
              if(metal.and.i_nu==0) then
                call phi_rpa_c_q0_met_sclp(nne,work,p00,p01,nrow,ncol,
     &                                     refrow,refcol)
              else
                call phi_rpa_c_q0_sclp(nne,i_nu,work,p_head(1,1,ind_nu),
     &                                 p_wing(1,1,ind_nu),nrow,ncol,
     &                                 refrow,refcol)
              endif
            endif
#endif
          else
            if(ncmpl==1) then
              if(metal.and.i_nu==0) then
                call phi_rpa_c_q0_met_r(n_pbt,nne,workr,p00,p01r)
              else
                call phi_rpa_c_q0_r(n_pbt,nne,i_nu,workr,
     &                              p_headr(1,1,ind_nu),
     &                              p_wingr(1,1,ind_nu))
              endif
            else
              if(metal.and.i_nu==0) then
                call phi_rpa_c_q0_met(n_pbt,nne,work,p00,p01)
              else
                call phi_rpa_c_q0(n_pbt,nne,i_nu,work,
     &                            p_head(1,1,ind_nu),p_wing(1,1,ind_nu))
              endif
            endif
          endif
        endif
c -------- Transform W to the original basis set -------------------
        if(lsclp) then
#ifdef SCLP
          if(ncmpl==1) then
            deallocate(workr)
            allocate(workr(nrow_pw(1),ncol_pw(1)))
            call pdgemm('t','n',n_pbt,nne,nne,1.d0,bm_sclpr,1,1,
     &                  desc_bm,w_sr,1,1,desc,0.d0,tmpr,1,1,desc_t)
            call pdgemm('n','n',n_pbt,n_pbt,nne,1.d0,tmpr,1,1,desc_t,
     &                  bm_sclpr,1,1,desc_bm,0.d0,workr,1,1,
     &                  desc_pw(1,1))
          else
            deallocate(work)
            allocate(work(nrow_pw(1),ncol_pw(1)))
            call pzgemm('c','n',n_pbt,nne,nne,(1.d0,0.d0),bm_sclp,1,1,
     &                  desc_bm,w_s,1,1,desc,(0.d0,0.d0),tmp,1,1,desc_t)
            call pzgemm('n','n',n_pbt,n_pbt,nne,(1.d0,0.d0),tmp,1,1,
     &                  desc_t,bm_sclp,1,1,desc_bm,(0.d0,0.d0),work,1,1,
     &                  desc_pw(1,1))
          endif
#endif
        else
          if(ncmpl==1) then
            call dgemm('t','n',n_pbt,nne,nne,1.d0,bm_coef_r,nne,
     &                 w_sr,nne,0.d0,tmpr,n_pbt)
            call dgemm('n','n',n_pbt,n_pbt,nne,1.d0,tmpr,n_pbt,
     &                 bm_coef_r,nne,0.d0,workr,n_pbt)
          else
            call zgemm('c','n',n_pbt,nne,nne,(1.d0,0.d0),bm_coef,nne,
     &                 w_s,nne,(0.d0,0.d0),tmp,n_pbt)
            call zgemm('n','n',n_pbt,n_pbt,nne,(1.d0,0.d0),tmp,n_pbt,
     &                 bm_coef,nne,(0.d0,0.d0),work,n_pbt)
          endif
        endif
c ------------------------------------------------------------------
        if(lsclp) then
          if(ncmpl==1) then
            call dlacpy('a',nrow_pw(1),ncol_pw(1),workr,nrow_pw(1),
     &                  pw_sclp(1,1,ind_nu,1),nrow_pw_max)
          else
            call zlacpy('a',nrow_pw(1),ncol_pw(1),work,nrow_pw(1),
     &                  pw_sclp(1,1,ind_nu,1),nrow_pw_max)
          endif
        else
          if(ncmpl==1) then
            call pack_hermit3r(workr,p_w_q_tau_mm(1,ind_nu,1),
     &                         p_w_q_tau_mi(1,1,ind_nu,1),
     &                         p_w_q_tau_ii(1,ind_nu,1),nplwgw(1),1)
          else
            call pack_hermit3(work,p_w_q_tau_mm(1,ind_nu,1),
     &                        p_w_q_tau_mi(1,1,ind_nu,1),
     &                        p_w_q_tau_ii(1,ind_nu,1),nplwgw(1))
          endif
        endif
c ------------------------------------------------------------------
        w_output(i_nu,1)=ws_head_nu(i_nu)*q2aver
        if(w_red) then
          if(lsclp) then
          else
            if(ncmpl==1) then
              allocate(work(n_pbt,n_pbt))
              call pb_c_from_r(n_pbt,n_pbt,n_pbt,work,workr,pnt(1,1))
            endif
            allocate(w_b(n_pbtot_red,n_pbtot_red))
            call rpb_from_fpb_full(w_b,work,1,n_pbtot_red,n_pbt)
            if(ncmpl==1) deallocate(work)
            n=n_pbmt_red+nplwgw_red(1)
            call pack_hermit(w_b,w_nu_q_b(1,1,ind_nu,1),n,
     &                       n_pbtot_red,n_pbtot_red,0.d0,1.d0)
            deallocate(w_b)
            allocate(w_b(nne,n_opt_pb_red))
            if(ncmpl==1) then
              allocate(w_s(nne,nne))
              w_s=w_sr
            endif
            call zgemm('n','n',nne,n_opt_pb_red,nne,(1.d0,0.d0),w_s,nne,
     &                 bb_red,nne,(0.d0,0.d0),w_b,nne)
            if(ncmpl==1) deallocate(w_s)
            call zgemm('c','n',n_opt_pb_red,n_opt_pb_red,nne,
     &                 (1.d0,0.d0),bb_red,nne,w_b,nne,(0.d0,0.d0),
     &                 w11_red(1,1,ind_nu),n_opt_pb_red)
            deallocate(w_b)
          endif
        endif
        if(ncmpl==1) deallocate(workr)
        if(ncmpl==2) deallocate(work)
c ------------------------------------------------------------------
      enddo  !! over ind_nu
      if(ncmpl==1) then
        deallocate(tmpr,w_sr,p_headr,p_wingr)
        if(metal) deallocate(p01r)
      else
        deallocate(tmp,w_s,p_head,p_wing)
        if(metal) deallocate(p01)
      endif
#ifdef SCLP
      if(lsclp) then
        if(ncmpl==1) deallocate(bm_sclpr)
        if(ncmpl==2) deallocate(bm_sclp)
      endif
#endif
      if(nproc_t/=1) call dgop(ws_head_nu,n_nu+1,'  +',comm_t)
      call nu_to_tau_spl(ws_head_nu,ws_head,1)
      ws_head=ws_head/2.d0
c   ------- Transform NU ---> TAU representation -------------------
      if(lsclp) then
        call nu_to_tau_ba3ut(pw_sclp,ncmpl*nrow_pw_max*ncol_pw_max,
     &                       ncmpl*nrow_pw_max*ncol_pw_max)
      else
        call nu_to_tau_ba3ut(p_w_q_tau_mm,ncmpl*ndim_pbmt(me_b+1),
     &                       ncmpl*ndim_pbmt(me_b+1))
        call nu_to_tau_ba3(p_w_q_tau_mi,ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                     ncmpl*nd_b_pbmt(me_b+1),nplwgw(1))
        call nu_to_tau_ba3ut(p_w_q_tau_ii,nn2,nn2)
      endif
      deallocate(vsq)
      end
