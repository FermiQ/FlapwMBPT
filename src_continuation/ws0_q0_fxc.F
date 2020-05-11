      subroutine ws0_q0_fxc(w_out_dir,e_out_dir,e1_out_dir)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(out) :: w_out_dir(n_phead_dir,ndim3_nu),
     &                       e_out_dir(n_phead_dir,ndim3_nu),
     &                       e1_out_dir(n_phead_dir,ndim3_nu)
      integer :: n_pbt,i,ind_nu,i_nu,nne,ii,j,nn2
      real*8 :: pi2a,pi8,p0,sqpi4,p00,r0
      complex*16 :: cc
      real*8, allocatable :: workr(:,:),w_sr(:,:),tmpr(:,:),
     &                       p_headr(:,:,:),p_wingr(:,:,:),
     &                       fxcr(:,:),xpr(:,:),fxc0r(:,:),f00(:),
     &                       x11(:)
      complex*16, allocatable :: w_s(:,:),tmp(:,:),work(:,:),
     &                           p_head(:,:,:),p_wing(:,:,:),
     &                           fxc(:,:),xp(:,:),fxc0(:,:)
      pi2a=(pi+pi)/par
      pi8=8.d0*pi
      sqpi4=sqrt(4.d0*pi)
      nn2=ncmpl*nd_b_pbint(me_b+1)
      n_pbt=n_pbmt+nplwgw(1)
      nne=n_opt_pb
c   ------- Transform P_FULL from TAU to NU representation for Q=0 -
      call tau_to_nu_ba3ut(p_w_q_tau_mm,ncmpl*ndim_pbmt(me_b+1),
     &                     ncmpl*ndim_pbmt(me_b+1),1)
      call tau_to_nu_ba3(p_w_q_tau_mi,ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                   ncmpl*nd_b_pbmt(me_b+1),nplwgw(1),1)
      call tau_to_nu_ba3ut(p_w_q_tau_ii,nn2,nn2,1)
c ------------------------------------------------------------------
      if(ncmpl==1) then
        allocate(w_sr(nne,nne))
        allocate(tmpr(n_pbt,nne))
      else
        allocate(w_s(nne,nne))
        allocate(tmp(n_pbt,nne))
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
        call output_p_tensor_r(p_headr)
        r0=0.d0
        if(me_t==0) then
          do i=1,3
            r0=r0+p_headr(i,i,1)
          enddo
          r0=r0/3.d0
        endif
        if(nproc_t/=1) call dgop(r0,1,'  +',comm_t)
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
        call output_p_tensor(p_head)
        r0=0.d0
        if(me_t==0) then
          do i=1,3
            r0=r0+p_head(i,i,1)
          enddo
          r0=r0/3.d0
        endif
        if(nproc_t/=1) call dgop(r0,1,'  +',comm_t)
      endif
c ---------- Get F_xc ----------------------------------------------
      allocate(x11(n_phead_dir))
      allocate(f00(n_phead_dir))
      do i=1,n_phead_dir
        f00(i)=pi8+fxc_00_bts(i)
      enddo
      if(ncmpl==1) then
        allocate(fxc0r(nne,nne))
        call unpack_fxc_r(nne,fxc0r,fxc_bts(1,1))
        do i=1,nne
          fxc0r(i,i)=v_opt_e(i)+fxc0r(i,i)
        enddo
        allocate(fxcr(nne,nne))
        fxcr=0.d0
      else
        allocate(fxc0(nne,nne))
        call unpack_fxc(nne,fxc0,fxc_bts(1,1))
        do i=1,nne
          fxc0(i,i)=v_opt_e(i)+fxc0(i,i)
        enddo
        allocate(fxc(nne,nne))
        fxc=(0.d0,0.d0)
      endif
c ------------------------------------------------------------------
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(ncmpl==1) then
          allocate(workr(n_pbt,n_pbt))
          call unpack_hermit3r(workr,p_w_q_tau_mm(1,ind_nu,1),
     &                         p_w_q_tau_mi(1,1,ind_nu,1),
     &                         p_w_q_tau_ii(1,ind_nu,1),nplwgw(1))
        else
          allocate(work(n_pbt,n_pbt))
          call unpack_hermit3(work,p_w_q_tau_mm(1,ind_nu,1),
     &                        p_w_q_tau_mi(1,1,ind_nu,1),
     &                        p_w_q_tau_ii(1,ind_nu,1),nplwgw(1))
        endif
        if(ncmpl==1) then
          call pw_from_pb_r(1,n_pbt,n_pbt,workr,p0,3)
        else
          call pw_from_pb(1,n_pbt,n_pbt,work,cc,3)
          p0=cc
        endif
c -------- Transform P to the orthonormal basis set ----------------
        if(ncmpl==1) then
          call dgemm('n','t',n_pbt,nne,n_pbt,1.d0,workr,n_pbt,
     &               bm_coef_r,nne,0.d0,tmpr,n_pbt)
          call dgemm('n','n',nne,nne,n_pbt,1.d0,bm_coef_r,nne,
     &               tmpr,n_pbt,0.d0,w_sr,nne)
          deallocate(workr)
        else
          call zgemm('n','c',n_pbt,nne,n_pbt,(1.d0,0.d0),work,n_pbt,
     &               bm_coef,nne,(0.d0,0.d0),tmp,n_pbt)
          call zgemm('n','n',nne,nne,n_pbt,(1.d0,0.d0),bm_coef,nne,
     &               tmp,n_pbt,(0.d0,0.d0),w_s,nne)
          deallocate(work)
        endif
        p_output(i_nu,1)=p0
c ---------- Evaluate R0 * (V + f) ---------------------------------
        if(ncmpl==1) then
          allocate(workr(nne,nne))
          call dgemm('n','n',nne,nne,nne,1.d0,w_sr,nne,fxc0r,nne,0.d0,
     &               workr,nne)
        else
          allocate(work(nne,nne))
          call zgemm('n','n',nne,nne,nne,(1.d0,0.d0),w_s,nne,fxc0,nne,
     &               (0.d0,0.d0),work,nne)
        endif
c ------------------------------------------------------------------
        do i=1,n_phead_dir
          x11(i)=1.d0-ph_dir(i,ind_nu)*f00(i)
        enddo
        if(ncmpl==1) then
          workr=-workr
          do i=1,nne
            workr(i,i)=1.d0+workr(i,i)
          enddo
          call invers_r(nne,workr,nne)
          allocate(xpr(nne,nne))
          call dgemm('n','n',nne,nne,nne,1.d0,workr,nne,w_sr,nne,0.d0,
     &               xpr,nne)
          deallocate(workr)
          do j=1,nne
            do i=1,nne
              w_sr(i,j)=v_opt_e(i)*xpr(i,j)*v_opt_e(j)
            enddo
          enddo
          deallocate(xpr)
        else
          work=-work
          do i=1,nne
            work(i,i)=(1.d0,0.d0)+work(i,i)
          enddo
          call invers_z(nne,work,nne)
          allocate(xp(nne,nne))
          call zgemm('n','n',nne,nne,nne,(1.d0,0.d0),work,nne,w_s,nne,
     &               (0.d0,0.d0),xp,nne)
          deallocate(work)
          do j=1,nne
            do i=1,nne
              w_s(i,j)=v_opt_e(i)*xp(i,j)*v_opt_e(j)
            enddo
          enddo
          deallocate(xp)
        endif
        do i=1,n_phead_dir
          p00=ph_dir(i,ind_nu)/x11(i)
          e1_out_dir(i,ind_nu)=1.d0+pi8*p00
          e_out_dir(i,ind_nu)=1.d0/e1_out_dir(i,ind_nu)
          w_out_dir(i,ind_nu)=e1_out_dir(i,ind_nu)*pi8
        enddo
        p00=0.d0
        do i=1,n_phead_dir
          p00=p00+w_head_dir(i)*ph_dir(i,ind_nu)/x11(i)
        enddo
        ws_head_nu(i_nu)=p00*pi8**2
        e_output(i_nu,1)=1.d0+pi8*p00
        w_output(i_nu,1)=ws_head_nu(i_nu)*q2aver
c -------- Transform W to the original basis set -------------------
        if(ncmpl==1) then
          call dgemm('t','n',n_pbt,nne,nne,1.d0,bm_coef_r,nne,
     &               w_sr,nne,0.d0,tmpr,n_pbt)
          allocate(workr(n_pbt,n_pbt))
          call dgemm('n','n',n_pbt,n_pbt,nne,1.d0,tmpr,n_pbt,
     &               bm_coef_r,nne,0.d0,workr,n_pbt)
        else
          call zgemm('c','n',n_pbt,nne,nne,(1.d0,0.d0),bm_coef,nne,
     &               w_s,nne,(0.d0,0.d0),tmp,n_pbt)
          allocate(work(n_pbt,n_pbt))
          call zgemm('n','n',n_pbt,n_pbt,nne,(1.d0,0.d0),tmp,n_pbt,
     &               bm_coef,nne,(0.d0,0.d0),work,n_pbt)
        endif
c ------------------------------------------------------------------
        if(ncmpl==1) then
          call pack_hermit3r(workr,p_w_q_tau_mm(1,ind_nu,1),
     &                       p_w_q_tau_mi(1,1,ind_nu,1),
     &                       p_w_q_tau_ii(1,ind_nu,1),nplwgw(1),1)
        else
          call pack_hermit3(work,p_w_q_tau_mm(1,ind_nu,1),
     &                      p_w_q_tau_mi(1,1,ind_nu,1),
     &                      p_w_q_tau_ii(1,ind_nu,1),nplwgw(1))
        endif
c ------------------------------------------------------------------
        if(ncmpl==1) deallocate(workr)
        if(ncmpl==2) deallocate(work)
c ------------------------------------------------------------------
        if(i_nu==0) then
          p0=0.d0
          do i=1,n_phead_dir
            p0=p0+w_head_dir(i)*ph_dir(i,ind_nu)
          enddo
          p0=pi8*p0
          fxc_00_bts=w_out_dir(:,1)/p0
          if(ncmpl==1) then
            do j=1,nne
              do i=1,nne
                p00=w_sr(i,j)
                if(i==j) p00=p00+v_opt_e(i)
                fxcr(i,j)=p00/p0
              enddo
            enddo
          else
            do j=1,nne
              do i=1,nne
                cc=w_s(i,j)
                if(i==j) cc=cc+v_opt_e(i)
                fxc(i,j)=cc/p0
              enddo
            enddo
          endif
        endif
      enddo  !! over ind_nu
      if(nproc_t*nproc_b/=1) then
        call brdcst(msgdbl,fxc_00_bts,8*n_phead_dir,0,comm_t_b)
        if(ncmpl==1) then
          call brdcst(msgdbl,fxcr,8*nne**2,0,comm_t_b)
        else
          call brdcst(msgdbl,fxc,16*nne**2,0,comm_t_b)
        endif
      endif
c ------------------ Mixing ----------------------------------------
      fxc_00_bts=0.3d0*fxc_00_bts+0.7d0*(f00-pi8)
      if(ncmpl==1) then
        do j=1,nne
          do i=1,nne
            r0=fxc0r(i,j)
            if(i==j) r0=r0-v_opt_e(i)
            fxcr(i,j)=0.3d0*fxcr(i,j)+0.7d0*r0
          enddo
        enddo
      else
        do j=1,nne
          do i=1,nne
            cc=fxc0(i,j)
            if(i==j) cc=cc-v_opt_e(i)
            fxc(i,j)=0.3d0*fxc(i,j)+0.7d0*cc
          enddo
        enddo
      endif
c ------------------------------------------------------------------
      if(ncmpl==1) then
        call pack_fxc_r(nne,fxcr,fxc_bts(1,1))
      else
        call pack_fxc(nne,fxc,fxc_bts(1,1))
      endif
      if(ncmpl==1) then
        deallocate(tmpr,w_sr,p_headr,p_wingr,fxcr,fxc0r)
      else
        deallocate(tmp,w_s,p_head,p_wing,fxc,fxc0)
      endif
      if(nproc_t/=1) call dgop(ws_head_nu,n_nu+1,'  +',comm_t)
      call nu_to_tau_spl(ws_head_nu,ws_head,1)
      ws_head=ws_head/2.d0
c   ------- Transform NU ---> TAU representation -------------------
      call nu_to_tau_ba3ut(p_w_q_tau_mm,ncmpl*ndim_pbmt(me_b+1),
     &                     ncmpl*ndim_pbmt(me_b+1))
      call nu_to_tau_ba3(p_w_q_tau_mi,ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                   ncmpl*nd_b_pbmt(me_b+1),nplwgw(1))
      call nu_to_tau_ba3ut(p_w_q_tau_ii,nn2,nn2)
      deallocate(f00,x11)
      end
