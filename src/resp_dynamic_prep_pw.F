      subroutine resp_dynamic_prep_pw(pw,key)
!     key = 1 for 00 Polarizability
!     key = 2 for -+ Polarizability
!     key = 3 for ZZ Polarizability
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      complex*16, intent(out) :: pw(0:n_nu,n_line)
      integer :: ind_k,ind_tau,k,iq,n_pbt,nn,ind_nu,i_nu,nn2
      complex*16 :: zdotc
      real*8, allocatable :: resp_0_tau_ii(:,:,:,:),tm1r(:),tmr(:,:),
     &                       resp_0_tau_mi(:,:,:,:,:),
     &                       resp_0_tau_mm(:,:,:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),vv(:)
      nn2=ncmpl*nd_b_pbint(me_b+1)
c --------- Response_0 on TAU ---------------------------------------
      allocate(resp_0_tau_mm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                       ndim3_k(me_k+1),nspin))
      allocate(resp_0_tau_mi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,ndim3_tn,
     &                       ndim3_k(me_k+1),nspin))
      allocate(resp_0_tau_ii(nn2,ndim3_tn,ndim3_k(me_k+1),nspin))
      resp_0_tau_mm=0.d0
      resp_0_tau_mi=0.d0
      resp_0_tau_ii=0.d0
      if(nspin==1) then
! in this case KEY doesn't have an effect for non-interacting case --      
        call pol_0_spin_spin(1,1,resp_0_tau_mm,resp_0_tau_mi,
     &                       resp_0_tau_ii)
        resp_0_tau_mm=resp_0_tau_mm*4.d0
        resp_0_tau_mi=resp_0_tau_mi*4.d0
        resp_0_tau_ii=resp_0_tau_ii*4.d0
      else if(nspin==2) then
        if(key==1.or.key==3) then
c ------- A(upup+dndn) and A(upup-dndn) -----------------------------
          call pol_0_spin_spin(2,2,resp_0_tau_mm(1,1,1,2),
     &                         resp_0_tau_mi(1,1,1,1,2),
     &                         resp_0_tau_ii(1,1,1,2)) ! up_up
          call pol_0_spin_spin(1,1,resp_0_tau_mm(1,1,1,1),
     &                         resp_0_tau_mi(1,1,1,1,1),
     &                         resp_0_tau_ii(1,1,1,1)) ! dn_dn
c ------------- Form proper combinations --------------------------
          do ind_k=1,ndim3_k(me_k+1)
            do ind_tau=1,ndim3_tau 
              resp_0_tau_mm(:,ind_tau,ind_k,1)=
     &	        (resp_0_tau_mm(:,ind_tau,ind_k,1)
     &	        +resp_0_tau_mm(:,ind_tau,ind_k,2))*2.d0 
              resp_0_tau_mi(:,:,ind_tau,ind_k,1)=
     &	        (resp_0_tau_mi(:,:,ind_tau,ind_k,1)
     &	        +resp_0_tau_mi(:,:,ind_tau,ind_k,2))*2.d0 
              resp_0_tau_ii(:,ind_tau,ind_k,1)=
     &	        (resp_0_tau_ii(:,ind_tau,ind_k,1)
     &	        +resp_0_tau_ii(:,ind_tau,ind_k,2))*2.d0
            enddo
          enddo
        else if(key==2) then
c ------- A(dnup+updn) and A(dnup-updn) ----------------------------
          call pol_0_spin_spin(2,1,resp_0_tau_mm(1,1,1,1),
     &                         resp_0_tau_mi(1,1,1,1,1),
     &                         resp_0_tau_ii(1,1,1,1)) ! up_dn
          call pol_0_spin_spin(1,2,resp_0_tau_mm(1,1,1,2),
     &                         resp_0_tau_mi(1,1,1,1,2),
     &                         resp_0_tau_ii(1,1,1,2)) ! dn_up
c ------------- Form proper A and B combinations for (+-) component
          allocate(tm1r(ncmpl*ndim_pbmt(me_b+1)))
          do ind_k=1,ndim3_k(me_k+1)
            do ind_tau=1,ndim3_tau
              tm1r=resp_0_tau_mm(:,ind_tau,ind_k,1)
     &	         +resp_0_tau_mm(:,ind_tau,ind_k,2)
              resp_0_tau_mm(:,ind_tau,ind_k,2)=
     &	          (resp_0_tau_mm(:,ind_tau,ind_k,1)
     &	          -resp_0_tau_mm(:,ind_tau,ind_k,2))*4.d0
              resp_0_tau_mm(:,ind_tau,ind_k,1)=tm1r*4.d0
            enddo
          enddo
          deallocate(tm1r)
          allocate(tmr(ncmpl*nd_b_pbmt(me_b+1),nplw_gw))
          do ind_k=1,ndim3_k(me_k+1)
            do ind_tau=1,ndim3_tau
              tmr=resp_0_tau_mi(:,:,ind_tau,ind_k,1)
     &	      +resp_0_tau_mi(:,:,ind_tau,ind_k,2)
              resp_0_tau_mi(:,:,ind_tau,ind_k,2)=
     &	          (resp_0_tau_mi(:,:,ind_tau,ind_k,1)
     &	          -resp_0_tau_mi(:,:,ind_tau,ind_k,2))*4.d0
              resp_0_tau_mi(:,:,ind_tau,ind_k,1)=tmr*4.d0
            enddo
          enddo
          deallocate(tmr)
          allocate(tm1r(nn2))
          do ind_k=1,ndim3_k(me_k+1)
            do ind_tau=1,ndim3_tau
              tm1r=resp_0_tau_ii(:,ind_tau,ind_k,1)
     &	          +resp_0_tau_ii(:,ind_tau,ind_k,2)
              resp_0_tau_ii(:,ind_tau,ind_k,2)=
     &	          (resp_0_tau_ii(:,ind_tau,ind_k,1)
     &	          -resp_0_tau_ii(:,ind_tau,ind_k,2))*4.d0
              resp_0_tau_ii(:,ind_tau,ind_k,1)=tm1r*4.d0
            enddo
          enddo
          deallocate(tm1r)
        endif
      endif
      pw=(0.d0,0.d0)
      nn=1
      if(key==2.and.nspin==2) nn=2
      allocate(vv(n_pbtot))
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        iq=k_npnt_in_line(k)
        if(iq==0) cycle
        n_pbt=n_pbmt+nplwgw(k)
        allocate(tmp(n_pbt,n_pbt))
        if(key==2) allocate(tmp1(n_pbt,n_pbt))
        call tau_to_nu_ba3ut(resp_0_tau_mm(1,1,ind_k,1),
     &	                     ncmpl*ndim_pbmt(me_b+1),
     &                       ncmpl*ndim_pbmt(me_b+1),1)
        call tau_to_nu_ba3(resp_0_tau_mi(1,1,1,ind_k,1),
     &	                   ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                     ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),1)
        call tau_to_nu_ba3ut(resp_0_tau_ii(1,1,ind_k,1),nn2,nn2,1)
        if(key==2.and.nspin==2) then
          call tau_to_nu_ba3ut(resp_0_tau_mm(1,1,ind_k,2),
     &                         ncmpl*ndim_pbmt(me_b+1),
     &                         ncmpl*ndim_pbmt(me_b+1),2)
          call tau_to_nu_ba3(resp_0_tau_mi(1,1,1,ind_k,2),
     &	                     ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                       ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),2)
          call tau_to_nu_ba3ut(resp_0_tau_ii(1,1,ind_k,2),nn2,nn2,2)
        endif
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          if(ncmpl==1) then
            call unpack_hermit3rOLD(tmp,resp_0_tau_mm(1,ind_nu,ind_k,1),
     &                              resp_0_tau_mi(1,1,ind_nu,ind_k,1),
     &                              resp_0_tau_ii(1,ind_nu,ind_k,1),
     &                              nplwgw(k),k)
          else
            call unpack_hermit3(tmp,resp_0_tau_mm(1,ind_nu,ind_k,1),
     &                          resp_0_tau_mi(1,1,ind_nu,ind_k,1),
     &                          resp_0_tau_ii(1,ind_nu,ind_k,1),
     &                          nplwgw(k))
          endif
          if(key==2.and.nspin==2) then
            if(ncmpl==1) then
              call unpack_hermit3rOLD(tmp1,
     &                                resp_0_tau_mm(1,ind_nu,ind_k,2),
     &                                resp_0_tau_mi(1,1,ind_nu,ind_k,2),
     &                                resp_0_tau_ii(1,ind_nu,ind_k,2),
     &                                nplwgw(k),k)
            else
              call unpack_hermit3(tmp1,resp_0_tau_mm(1,ind_nu,ind_k,2),
     &                            resp_0_tau_mi(1,1,ind_nu,ind_k,2),
     &                            resp_0_tau_ii(1,ind_nu,ind_k,2),
     &                            nplwgw(k))
            endif
            tmp=tmp+(0.d0,1.d0)*tmp1
          endif
          call zgemv('n',n_pbt,n_pbt,(1.d0,0.d0),tmp,n_pbt,
     &	             pw_pb(1,ind_k),1,(0.d0,0.d0),vv,1)
          pw(i_nu,iq)=zdotc(n_pbt,pw_pb(1,ind_k),1,vv,1)
        enddo
        deallocate(tmp)
        if(key==2) deallocate(tmp1)
      enddo
      deallocate(resp_0_tau_mm,resp_0_tau_mi,resp_0_tau_ii,vv)
      if(nproc_t/=1) call dgop(pw,2*(n_nu+1)*n_line,'  +',comm_t)
      if(nproc_k/=1) call dgop(pw,2*(n_nu+1)*n_line,'  +',comm_k)
      call timel('** RESP_DYNAMIC_PREP_PW finished ***')
      end
