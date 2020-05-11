#ifdef SCLP
      subroutine p_full_lda_sclp(ch,ind_k,n_ch0,ind_ch0,n_ix0,ind_ix0)
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
      character*5, intent(in) :: ch
      integer, intent(in) :: ind_k,n_ch0,ind_ch0(n_ch0),n_ix0,
     &                       ind_ix0(n_ix0)
      integer :: k,n_pbt,i,ind_nu,i_nu,info,n,i00,i0z,izz,nrow,ncol,j,
     &           jj,ii
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: tmpr(:,:),thet_q(:,:,:),r00r(:,:),
     &                       r0zr(:,:),a11r(:,:),a12r(:,:),a21r(:,:),
     &                       a22r(:,:)
      complex*16, allocatable :: tmp(:,:),r00(:,:),r0z(:,:),a11(:,:),
     &                           a12(:,:),a21(:,:),a22(:,:)
      k=n3_mpi_k(me_k+1)+ind_k
      n_pbt=n_pbmt+nplwgw(k)
      nrow=nrow_pw(k)
      ncol=ncol_pw(k)
      if(ch=='___00') n=1
      if(ch=='00_ZZ') n=2
      allocate(ipiv(nrow+nblock_pw))
      i00=ind_ixc(1)
      if(n==2) then
        i0z=ind_ixc(3)
        izz=ind_ixc(2)
      endif
c ------------------------------------------------------------------
      allocate(thet_q(ncmpl*nrow,ncol,n_ixc_0))
      call theta_q_ss_sclp(k,nrow,ncol,thet_q)
      if(k==1) call timel('**** F-xc finished *****************')
c   --- Transform P_FULL from TAU to NU representation for given Q -
      call tau_to_nu_ba3ut(pw_sclp(1,1,1,ind_k),
     &                     ncmpl*nrow_pw_max*ncol_pw_max,
     &                     ncmpl*nrow_pw_max*ncol_pw_max,1)
      if(n==2) then
        call tau_to_nu_ba3ut(pz_sclp(1,1,1,ind_k),
     &                       ncmpl*nrow_pw_max*ncol_pw_max,
     &                       ncmpl*nrow_pw_max*ncol_pw_max,1)
      endif
c ------------------------------------------------------------------
      if(k==1) call timel('**** Nu-loop starts ****************')
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(ncmpl==1) then
          allocate(r00r(nrow,ncol))
          call dlacpy('a',nrow,ncol,pw_sclp(1,1,ind_nu,ind_k),
     &                nrow_pw_max,r00r,nrow)
          allocate(a11r(nrow,ncol))
          call pdgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r00r,1,1,
     &                desc_pw(1,k),thet_q,1,1,desc_pw(1,k),0.d0,a11r,1,
     &                1,desc_pw(1,k))
          do j=1,ncol
            jj=refcol_pw(j,k)
            do i=1,nrow
              ii=refrow_pw(i,k)
              if(ii==jj) a11r(i,j)=1.d0+a11r(i,j)
            enddo
          enddo
          if(n==1) then
            call pdgesv(n_pbt,n_pbt,a11r,1,1,desc_pw(1,k),ipiv,r00r,1,1,
     &                  desc_pw(1,k),info)
            deallocate(a11r)
          else
            allocate(r0zr(nrow,ncol))
            call dlacpy('a',nrow,ncol,pz_sclp(1,1,ind_nu,ind_k),
     &                  nrow_pw_max,r0zr,nrow)
            allocate(a22r(nrow,ncol))
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r0zr,1,1,
     &                  desc_pw(1,k),thet_q(1,1,i0z),1,1,desc_pw(1,k),
     &                  0.d0,a22r,1,1,desc_pw(1,k))
            a11r=a11r+a22r
            do j=1,ncol
              jj=refcol_pw(j,k)
              do i=1,nrow
                ii=refrow_pw(i,k)
                if(ii==jj) a22r(i,j)=1.d0+a22r(i,j)
              enddo
            enddo
            allocate(a12r(nrow,ncol))
            allocate(a21r(nrow,ncol))
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r00r,1,1,
     &                  desc_pw(1,k),thet_q(1,1,i0z),1,1,desc_pw(1,k),
     &                  0.d0,a12r,1,1,desc_pw(1,k))
            a21r=a12r
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r0zr,1,1,
     &                  desc_pw(1,k),thet_q(1,1,izz),1,1,desc_pw(1,k),
     &                  1.d0,a12r,1,1,desc_pw(1,k))
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r0zr,1,1,
     &                  desc_pw(1,k),thet_q(1,1,i00),1,1,desc_pw(1,k),
     &                  1.d0,a21r,1,1,desc_pw(1,k))
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,r00r,1,1,
     &                  desc_pw(1,k),thet_q(1,1,izz),1,1,desc_pw(1,k),
     &                  1.d0,a22r,1,1,desc_pw(1,k))
            call invers_r_sclp(n_pbt,nrow,ncol,nblock_pw,desc_pw(1,k),
     &                         a22r)
            allocate(tmpr(nrow,ncol))
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,1.d0,a12r,1,1,
     &                  desc_pw(1,k),a22r,1,1,desc_pw(1,k),0.d0,tmpr,1,
     &                  1,desc_pw(1,k))
            deallocate(a12r,a22r)
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,tmpr,1,1,
     &                  desc_pw(1,k),a21r,1,1,desc_pw(1,k),1.d0,a11r,1,
     &                  1,desc_pw(1,k))
            call pdgemm('n','n',n_pbt,n_pbt,n_pbt,-1.d0,tmpr,1,1,
     &                  desc_pw(1,k),r0zr,1,1,desc_pw(1,k),1.d0,r00r,1,
     &                  1,desc_pw(1,k))
            deallocate(tmpr,a21r,r0zr)
            call pdgesv(n_pbt,n_pbt,a11r,1,1,desc_pw(1,k),ipiv,r00r,1,1,
     &                  desc_pw(1,k),info)
            deallocate(a11r)
          endif
          call dlacpy('a',nrow,ncol,r00r,nrow,pw_sclp(1,1,ind_nu,ind_k),
     &                nrow_pw_max)
          deallocate(r00r)
        else if(ncmpl==2) then
          allocate(r00(nrow,ncol))
          call zlacpy('a',nrow,ncol,pw_sclp(1,1,ind_nu,ind_k),
     &                nrow_pw_max,r00,nrow)
          allocate(a11(nrow,ncol))
          call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r00,1,1,
     &                desc_pw(1,k),thet_q,1,1,desc_pw(1,k),(0.d0,0.d0),
     &                a11,1,1,desc_pw(1,k))
          do j=1,ncol
            jj=refcol_pw(j,k)
            do i=1,nrow
              ii=refrow_pw(i,k)
              if(ii==jj) a11(i,j)=(1.d0,0.d0)+a11(i,j)
            enddo
          enddo
          if(n==1) then
            call pzgesv(n_pbt,n_pbt,a11,1,1,desc_pw(1,k),ipiv,r00,1,1,
     &                  desc_pw(1,k),info)
            deallocate(a11)
          else
            allocate(r0z(nrow,ncol))
            call zlacpy('a',nrow,ncol,pz_sclp(1,1,ind_nu,ind_k),
     &                  nrow_pw_max,r0z,nrow)
            allocate(a22(nrow,ncol))
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r0z,1,1,
     &                  desc_pw(1,k),thet_q(1,1,i0z),1,1,desc_pw(1,k),
     &                  (0.d0,0.d0),a22,1,1,desc_pw(1,k))
            a11=a11+a22
            do j=1,ncol
              jj=refcol_pw(j,k)
              do i=1,nrow
                ii=refrow_pw(i,k)
                if(ii==jj) a22(i,j)=(1.d0,0.d0)+a22(i,j)
              enddo
            enddo
            allocate(a12(nrow,ncol))
            allocate(a21(nrow,ncol))
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r00,1,1,
     &                  desc_pw(1,k),thet_q(1,1,i0z),1,1,desc_pw(1,k),
     &                  (0.d0,0.d0),a12,1,1,desc_pw(1,k))
            a21=a12
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r0z,1,1,
     &                  desc_pw(1,k),thet_q(1,1,izz),1,1,desc_pw(1,k),
     &                  (1.d0,0.d0),a12,1,1,desc_pw(1,k))
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r0z,1,1,
     &                  desc_pw(1,k),thet_q(1,1,i00),1,1,desc_pw(1,k),
     &                  (1.d0,0.d0),a21,1,1,desc_pw(1,k))
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),r00,1,1,
     &                  desc_pw(1,k),thet_q(1,1,izz),1,1,desc_pw(1,k),
     &                  (1.d0,0.d0),a22,1,1,desc_pw(1,k))
            call invers_z_sclp(n_pbt,nrow,ncol,nblock_pw,desc_pw(1,k),
     &                         a22)
            allocate(tmp(nrow,ncol))
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),a12,1,1,
     &                  desc_pw(1,k),a22,1,1,desc_pw(1,k),(0.d0,0.d0),
     &                  tmp,1,1,desc_pw(1,k))
            deallocate(a12,a22)
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),tmp,1,1,
     &                  desc_pw(1,k),a21,1,1,desc_pw(1,k),(1.d0,0.d0),
     &                  a11,1,1,desc_pw(1,k))
            call pzgemm('n','n',n_pbt,n_pbt,n_pbt,(-1.d0,0.d0),tmp,1,1,
     &                  desc_pw(1,k),r0z,1,1,desc_pw(1,k),(1.d0,0.d0),
     &                  r00,1,1,desc_pw(1,k))
            deallocate(tmp,a21,r0z)
            call pzgesv(n_pbt,n_pbt,a11,1,1,desc_pw(1,k),ipiv,r00,1,1,
     &                  desc_pw(1,k),info)
            deallocate(a11)
          endif
          call zlacpy('a',nrow,ncol,r00,nrow,pw_sclp(1,1,ind_nu,ind_k),
     &                nrow_pw_max)
          deallocate(r00)
        endif
      enddo  !! over ind_nu
      if(k==1) call timel('**** Nu-loop finished **************')
c   ------- Transform NU ---> TAU representation -------------------
      call nu_to_tau_ba3ut(pw_sclp(1,1,1,ind_k),
     &                     ncmpl*nrow_pw_max*ncol_pw_max,
     &                     ncmpl*nrow_pw_max*ncol_pw_max)
      deallocate(ipiv)
      if(k==1) call timel('**** P_FULL_LDA_SCLP finished ******')
      end
#endif
