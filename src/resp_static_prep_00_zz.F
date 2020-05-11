      subroutine resp_static_prep_00_zz(p_non)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(out) :: p_non(n_pbtot,n_pbtot,nc_line,nspin)
      integer :: ind_k,ind_tau,ispin,k,iq,n_pbt,kc,j1,i1,nn2
      real*8 :: tmp
      complex*16 :: tmpc
      complex*16, allocatable :: ppp(:,:)
      nn2=ncmpl*nd_b_pbint(me_b+1)
c --------- Response_0 on TAU ---------------------------------------
      allocate(resp_lda_tau_mm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                         ndim3_k(me_k+1),nspin))
      allocate(resp_lda_tau_mi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                         ndim3_tn,ndim3_k(me_k+1),nspin))
      allocate(resp_lda_tau_ii(nn2,ndim3_tn,ndim3_k(me_k+1),nspin))
      resp_lda_tau_mm=0.d0
      resp_lda_tau_mi=0.d0
      resp_lda_tau_ii=0.d0
      if(nspin==1) then
        call pol_0_spin_spin(1,1,resp_lda_tau_mm,resp_lda_tau_mi,
     &                       resp_lda_tau_ii)
        resp_lda_tau_mm=resp_lda_tau_mm*4.d0
        resp_lda_tau_mi=resp_lda_tau_mi*4.d0
        resp_lda_tau_ii=resp_lda_tau_ii*4.d0
      else if(nspin==2) then
c ------- A(upup+dndn) and A(upup-dndn) -----------------------------
        call pol_0_spin_spin(2,2,resp_lda_tau_mm(1,1,1,2),
     &                       resp_lda_tau_mi(1,1,1,1,2),
     &                       resp_lda_tau_ii(1,1,1,2)) ! up_up
        call pol_0_spin_spin(1,1,resp_lda_tau_mm(1,1,1,1),
     &                       resp_lda_tau_mi(1,1,1,1,1),
     &                       resp_lda_tau_ii(1,1,1,1)) ! dn_dn
c ------------- Form proper combinations --------------------------
        do ind_k=1,ndim3_k(me_k+1)
          do ind_tau=1,ndim3_tau
            do i1=1,ncmpl*ndim_pbmt(me_b+1)
              tmp=resp_lda_tau_mm(i1,ind_tau,ind_k,1)
     &	        +resp_lda_tau_mm(i1,ind_tau,ind_k,2)
              resp_lda_tau_mm(i1,ind_tau,ind_k,2)=
     &	        (resp_lda_tau_mm(i1,ind_tau,ind_k,2)
     &	        -resp_lda_tau_mm(i1,ind_tau,ind_k,1))*2.d0
              resp_lda_tau_mm(i1,ind_tau,ind_k,1)=tmp*2.d0
            enddo
            do j1=1,nplw_gw
              do i1=1,ncmpl*nd_b_pbmt(me_b+1)
                tmp=resp_lda_tau_mi(i1,j1,ind_tau,ind_k,1)
     &	           +resp_lda_tau_mi(i1,j1,ind_tau,ind_k,2)
                resp_lda_tau_mi(i1,j1,ind_tau,ind_k,2)=
     &	        (resp_lda_tau_mi(i1,j1,ind_tau,ind_k,2)
     &	        -resp_lda_tau_mi(i1,j1,ind_tau,ind_k,1))*2.d0
                resp_lda_tau_mi(i1,j1,ind_tau,ind_k,1)=tmp*2.d0
              enddo
            enddo
            do i1=1,nn2
              tmp=resp_lda_tau_ii(i1,ind_tau,ind_k,1)
     &	         +resp_lda_tau_ii(i1,ind_tau,ind_k,2)
              resp_lda_tau_ii(i1,ind_tau,ind_k,2)=
     &	        (resp_lda_tau_ii(i1,ind_tau,ind_k,2)
     &	        -resp_lda_tau_ii(i1,ind_tau,ind_k,1))*2.d0
              resp_lda_tau_ii(i1,ind_tau,ind_k,1)=tmp*2.d0
            enddo
          enddo
        enddo
      endif
      p_non=(0.d0,0.d0)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
          kc=k_c_from_a(k)
          if(kc==0) cycle
          iq=k_npnt_in_line_c(kc)
          if(iq==0) cycle
          n_pbt=n_pbmt+nplwgw(k)
          call tau_to_nu_ba3ut(resp_lda_tau_mm(1,1,ind_k,ispin),
     &	                       ncmpl*ndim_pbmt(me_b+1),
     &                         ncmpl*ndim_pbmt(me_b+1),1)
          call tau_to_nu_ba3(resp_lda_tau_mi(1,1,1,ind_k,ispin),
     &	                     ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                       ncmpl*nd_b_pbmt(me_b+1),nplwgw(k),1)
          call tau_to_nu_ba3ut(resp_lda_tau_ii(1,1,ind_k,ispin),nn2,
     &                         nn2,1)
          if(me_t==0) then
            allocate(ppp(n_pbt,n_pbt))
            if(ncmpl==1) then
              call unpack_hermit3rOLD(ppp,
     &                                resp_lda_tau_mm(1,1,ind_k,ispin),
     &                               resp_lda_tau_mi(1,1,1,ind_k,ispin),
     &                                resp_lda_tau_ii(1,1,ind_k,ispin),
     &                                nplwgw(k),k)
            else
              call unpack_hermit3(ppp,resp_lda_tau_mm(1,1,ind_k,ispin),
     &                            resp_lda_tau_mi(1,1,1,ind_k,ispin),
     &                            resp_lda_tau_ii(1,1,ind_k,ispin),
     &                            nplwgw(k))
            endif
            p_non(1:n_pbt,1:n_pbt,iq,ispin)=ppp
            deallocate(ppp)
          endif
        enddo
      enddo
      deallocate(resp_lda_tau_mm,resp_lda_tau_mi,resp_lda_tau_ii)
      if(nproc_t/=1) call dgop(p_non,2*n_pbtot**2*nc_line*nspin,'  +',
     &	                         comm_t)
      if(nproc_k/=1) call dgop(p_non,2*n_pbtot**2*nc_line*nspin,'  +',
     &	                       comm_k)
      call timel('** RESP_STATIC_PREP_00_ZZ finished *')
      end
