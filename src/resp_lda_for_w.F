      subroutine resp_lda_for_w
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
      integer :: ind_k,ind_tau,j1,i1,nn2,k
      real*8 :: tmp
      nn2=ncmpl*nd_b_pbint(me_b+1)
      if(nspin==1) then
        n_chi_0=1
        ind_chi(1)=1
        n_ixc_0=1
        ind_ixc(n_ixc_0)=1
      else if(nspin==2) then
        n_ixc_0=3
        ind_ixc(1)=1
        ind_ixc(2)=2
        ind_ixc(3)=3
        n_chi_0=2
        ind_chi(1)=1
        ind_chi(2)=2
      endif
c ------- Exchange-Correlation Kernel ------------------------------
      call theta_xc(0)
      if(lsclp) then
        pw_sclp=0.d0
      else
        p_w_q_tau_mm=0.d0
        p_w_q_tau_mi=0.d0
        p_w_q_tau_ii=0.d0
      endif
      if(nspin==2) then
        if(lsclp) then
          allocate(pz_sclp(ncmpl*nrow_pw_max,ncol_pw_max,ndim3_tn,
     &                     ndim3_k(me_k+1)))
          pz_sclp=0.d0
        else
          allocate(p_0z_mm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                     ndim3_k(me_k+1)))
          allocate(p_0z_mi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,ndim3_tn,
     &                     ndim3_k(me_k+1)))
          allocate(p_0z_ii(nn2,ndim3_tn,ndim3_k(me_k+1)))
          p_0z_mm=0.d0
          p_0z_mi=0.d0
          p_0z_ii=0.d0
        endif
      endif
      if(nspin==1) then
        if(lsclp) then
          call pol_sclp_spin_spin(1,1,pw_sclp)
          pw_sclp=pw_sclp*4.d0/nrel
        else
          call pol_0_spin_spin(1,1,p_w_q_tau_mm,p_w_q_tau_mi,
     &                         p_w_q_tau_ii)
          p_w_q_tau_mm=p_w_q_tau_mm*4.d0/nrel
          p_w_q_tau_mi=p_w_q_tau_mi*4.d0/nrel
          p_w_q_tau_ii=p_w_q_tau_ii*4.d0/nrel
        endif
      else if(nspin==2) then
        if(lsclp) then
c ------- A(upup+dndn) and A(upup-dndn) ----------------------------
          call pol_sclp_spin_spin(2,2,pw_sclp)    ! up_up
          call pol_sclp_spin_spin(1,1,pz_sclp)    ! dn_dn
c ------------- Form proper combinations ---------------------------
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            do ind_tau=1,ndim3_tau
              do j1=1,ncmpl*ncol_pw(k)
                do i1=1,nrow_pw(k)
                  tmp=pw_sclp(i1,j1,ind_tau,ind_k)
     &               +pz_sclp(i1,j1,ind_tau,ind_k)
                  pz_sclp(i1,j1,ind_tau,ind_k)=
     &                (pw_sclp(i1,j1,ind_tau,ind_k)
     &                -pz_sclp(i1,j1,ind_tau,ind_k))*2.d0
                  pw_sclp(i1,j1,ind_tau,ind_k)=tmp*2.d0
                enddo
              enddo
            enddo
          enddo
        else
c ------- A(upup+dndn) and A(upup-dndn) ----------------------------
          call pol_0_spin_spin(2,2,p_w_q_tau_mm,p_w_q_tau_mi,
     &                         p_w_q_tau_ii)                   ! up_up
          call pol_0_spin_spin(1,1,p_0z_mm,p_0z_mi,p_0z_ii)    ! dn_dn
c ------------- Form proper combinations ---------------------------
          do ind_k=1,ndim3_k(me_k+1)
            do ind_tau=1,ndim3_tau
              do i1=1,ncmpl*ndim_pbmt(me_b+1)
                tmp=p_w_q_tau_mm(i1,ind_tau,ind_k)
     &             +p_0z_mm(i1,ind_tau,ind_k)
                p_0z_mm(i1,ind_tau,ind_k)=
     &                (p_w_q_tau_mm(i1,ind_tau,ind_k)
     &                -p_0z_mm(i1,ind_tau,ind_k))*2.d0
                p_w_q_tau_mm(i1,ind_tau,ind_k)=tmp*2.d0
              enddo
              do j1=1,nplw_gw
                do i1=1,ncmpl*nd_b_pbmt(me_b+1)
                  tmp=p_w_q_tau_mi(i1,j1,ind_tau,ind_k)
     &               +p_0z_mi(i1,j1,ind_tau,ind_k)
                  p_0z_mi(i1,j1,ind_tau,ind_k)=
     &                (p_w_q_tau_mi(i1,j1,ind_tau,ind_k)
     &                -p_0z_mi(i1,j1,ind_tau,ind_k))*2.d0
                  p_w_q_tau_mi(i1,j1,ind_tau,ind_k)=tmp*2.d0
                enddo
              enddo
              do i1=1,nn2
                tmp=p_w_q_tau_ii(i1,ind_tau,ind_k)
     &             +p_0z_ii(i1,ind_tau,ind_k)
                p_0z_ii(i1,ind_tau,ind_k)=
     &                (p_w_q_tau_ii(i1,ind_tau,ind_k)
     &                -p_0z_ii(i1,ind_tau,ind_k))*2.d0
                p_w_q_tau_ii(i1,ind_tau,ind_k)=tmp*2.d0
              enddo
            enddo
          enddo
        endif
      endif
      call timel('**** RESP_LDA_FOR_W finished *******')
      end
