      subroutine p_q0_num_prep(ph,mm,mi,ii)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: ii(ncmpl*nd_b_pbint(me_b+1),ndim3_tn,
     &                         ndim3_k(me_k+1)),
     &                      mi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                         ndim3_tn,ndim3_k(me_k+1)),
     &                      mm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                         ndim3_k(me_k+1))
      real*8, intent(out) :: ph(n_k_head_0,ndim3_nu)
      integer :: n_pbt,i,i1,ind_q,iq,npw,ind_tau
      complex*16 :: zdotc
      real*8, allocatable :: pht(:,:),workr(:,:),tmpr(:)
      complex*16, allocatable :: work(:,:),tmp(:)
      allocate(pht(n_k_head_0,ndim3_tau))
      pht=0.d0
      do ind_q=1,ndim3_k(me_k+1)
        iq=n3_mpi_k(me_k+1)+ind_q
        i1=0
        do i=1,n_k_head_0
          if(iq==khead0(i)) i1=i
        enddo
        if(i1==0) cycle
        npw=nplwgw(iq)
        n_pbt=n_pbmt+npw
        if(ncmpl==1) allocate(workr(n_pbt,n_pbt),tmpr(n_pbt))
        if(ncmpl==2) allocate(work(n_pbt,n_pbt),tmp(n_pbt))
        do ind_tau=1,ndim3_tau
          if(ncmpl==1) then
            workr=0.d0
            call boson_unpack_tau3allr(workr,mm(1,ind_tau,ind_q),
     &                                 mi(1,1,ind_tau,ind_q),
     &                                 ii(1,ind_tau,ind_q),npw)
            call dgemv('n',n_pbt,n_pbt,1.d0,workr,n_pbt,
     &                 pw_pb_tild_r(1,ind_q),1,0.d0,tmpr,1)
            pht(i1,ind_tau)=dot_product(pw_pb_tild_r(1:n_pbt,ind_q),
     &                                  tmpr)
          else
            work=(0.d0,0.d0)
            call boson_unpack_tau3all(work,mm(1,ind_tau,ind_q),
     &                                mi(1,1,ind_tau,ind_q),
     &                                ii(1,ind_tau,ind_q),npw)
            call zgemv('n',n_pbt,n_pbt,(1.d0,0.d0),work,n_pbt,
     &                 pw_pb_tild(1,ind_q),1,(0.d0,0.d0),tmp,1)
            pht(i1,ind_tau)=zdotc(n_pbt,pw_pb_tild(1,ind_q),1,tmp,1)
          endif
        enddo
        if(ncmpl==1) deallocate(workr,tmpr)
        if(ncmpl==2) deallocate(work,tmp)
      enddo
      if(nproc_k/=1) then
        call dgop(pht,n_k_head_0*ndim3_tau,'  +',comm_k)
      endif
      call tau_to_nu_aa_g(ph,pht,n_k_head_0)
      deallocate(pht)
      end



      subroutine p_q0_num_prep_sclp(ph,pwsclp)
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
      real*8, intent(in) :: pwsclp(ncmpl*nrow_pw_max,ncol_pw_max,
     &                             ndim3_tn,ndim3_k(me_k+1))
      real*8, intent(out) :: ph(n_k_head_0,ndim3_nu)
      integer :: i,i1,ind_q,iq,ind_tau,nrow,ncol,ii,j,jj
      complex*16 :: zdotc
      real*8, allocatable :: pht(:,:),workr(:,:)
      complex*16, allocatable :: work(:,:)
      allocate(pht(n_k_head_0,ndim3_tau))
      pht=0.d0
      do ind_q=1,ndim3_k(me_k+1)
        iq=n3_mpi_k(me_k+1)+ind_q
        i1=0
        do i=1,n_k_head_0
          if(iq==khead0(i)) i1=i
        enddo
        if(i1==0) cycle
        nrow=nrow_pw(iq)
        ncol=ncol_pw(iq)
        if(ncmpl==1) allocate(workr(nrow,ncol))
        if(ncmpl==2) allocate(work(nrow,ncol))
        do ind_tau=1,ndim3_tau
          if(ncmpl==1) then
            call dlacpy('a',nrow,ncol,pwsclp(1,1,ind_tau,ind_q),
     &                  nrow_pw_max,workr,nrow)
            workr=workr/2
            do j=1,ncol
              jj=refcol_pw(j,iq)
              do i=1,nrow
                ii=refrow_pw(i,iq)
                pht(i1,ind_tau)=pht(i1,ind_tau)+pw_pb_tild_r(ii,ind_q)
     &                                         *workr(i,j)
     &                                         *pw_pb_tild_r(jj,ind_q)
              enddo
            enddo
          else
            call zlacpy('a',nrow,ncol,pwsclp(1,1,ind_tau,ind_q),
     &                  nrow_pw_max,work,nrow)
            work=work/2
            do j=1,ncol
              jj=refcol_pw(j,iq)
              do i=1,nrow
                ii=refrow_pw(i,iq)
                pht(i1,ind_tau)=pht(i1,ind_tau)
     &                         +conjg(pw_pb_tild(ii,ind_q))
     &                         *work(i,j)
     &                         *pw_pb_tild(jj,ind_q)
              enddo
            enddo
          endif
        enddo
        if(ncmpl==1) deallocate(workr)
        if(ncmpl==2) deallocate(work)
      enddo
      call dgop(pht,n_k_head_0*ndim3_tau,'  +',comm_b)
      if(nproc_k/=1) then
        call dgop(pht,n_k_head_0*ndim3_tau,'  +',comm_k)
      endif
      call tau_to_nu_aa_g(ph,pht,n_k_head_0)
      deallocate(pht)
      end
