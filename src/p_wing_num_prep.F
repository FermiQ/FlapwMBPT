      subroutine p_wing_num_prep(pwng,mm,mi,ii)
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
      real*8, intent(out) :: pwng(ncmpl*n_pbtot,0:3,ndim3_nu)
      integer :: n_pbt,i,i1,ind_q,iq,npw,ind_tau,iqr(0:3),kq
      real*8, allocatable :: workr(:,:),tmpr(:),pwt(:,:,:)
      complex*16, allocatable :: work(:,:),tmp(:)
      allocate(pwt(ncmpl*n_pbtot,0:3,ndim3_tau))
      pwt=0.d0
      iqr(0)=1
      do i1=1,3
        call zone1_number(qb0(1,i1),rb0,ndiv,kq)
        kq=index_k1(kq)
        iqr(i1)=i_kref(kq)
      enddo
c ------------------------------------------------------------------
      do ind_q=1,ndim3_k(me_k+1)
        iq=n3_mpi_k(me_k+1)+ind_q
        do i=0,3
          if(iq/=iqr(i)) cycle
          npw=nplwgw(iq)
          n_pbt=n_pbmt+npw
          if(ncmpl==1) allocate(workr(n_pbt,n_pbt),tmpr(n_pbt))
          if(ncmpl==2) allocate(work(n_pbt,n_pbt),tmp(n_pbt))
          do ind_tau=1,ndim3_tau
            if(ncmpl==1) then
              workr=0.d0
              call boson_unpack_tau3allr(workr,mm(1,ind_tau,ind_q),
     &                                   mi(1,1,ind_tau,ind_q),
     &                                   ii(1,ind_tau,ind_q),npw)
              call dgemv('n',n_pbt,n_pbt,1.d0,workr,n_pbt,
     &                   pw_pb_tild_r(1,ind_q),1,0.d0,tmpr,1)
              pwt(1:n_pbt,i,ind_tau)=tmpr
            else
              work=(0.d0,0.d0)
              call boson_unpack_tau3all(work,mm(1,ind_tau,ind_q),
     &                                  mi(1,1,ind_tau,ind_q),
     &                                  ii(1,ind_tau,ind_q),npw)
              call zgemv('n',n_pbt,n_pbt,(1.d0,0.d0),work,n_pbt,
     &                   pw_pb_tild(1,ind_q),1,(0.d0,0.d0),tmp,1)
              tmp=conjg(tmp)
              call zcopy(n_pbt,tmp,1,pwt(1,i,ind_tau),1)
            endif
          enddo
          if(ncmpl==1) deallocate(workr,tmpr)
          if(ncmpl==2) deallocate(work,tmp)
        enddo
      enddo
      if(nproc_k/=1) then
        call dgop(pwt,ncmpl*n_pbtot*4*ndim3_tau,'  +',comm_k)
      endif
      call tau_to_nu_aa_g(pwng,pwt,ncmpl*n_pbtot*4)
      deallocate(pwt)
      end
