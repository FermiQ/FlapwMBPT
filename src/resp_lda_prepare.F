      subroutine resp_lda_prepare(chi_tmp,cht)
c     cht = T if we need to store P as a function of TAU      
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical, intent(in) :: cht
      character*3, intent(in) :: chi_tmp
      integer :: ind_k,i,j,ind_tau,j1,i1,nn2
      real*8 :: tmp
      complex*16 :: tmpc
      nn2=ncmpl*nd_b_pbint(me_b+1)
      n_ixc_0=0
      n_chi_0=0
      ind_chi=0
      if(nspin==1) then
        n_chi_0=1
        ind_chi(1)=1
        if(chi_tmp(1:1)=='1') then
          n_ixc_0=n_ixc_0+1
          ind_ixc(n_ixc_0)=1
        endif
        if(chi_tmp(3:3)=='1') then
          n_ixc_0=n_ixc_0+1
          ind_ixc(n_ixc_0)=2
        endif
      else if(nspin==2) then
        if(chi_tmp(1:1)=='1'.or.chi_tmp(3:3)=='1') then
          n_ixc_0=n_ixc_0+3
          ind_ixc(n_ixc_0-2)=1
          ind_ixc(n_ixc_0-1)=2
          ind_ixc(n_ixc_0)=3
          n_chi_0=n_chi_0+2
          ind_chi(1)=1
          ind_chi(2)=2
        endif
        if(chi_tmp(2:2)=='1') then
          n_ixc_0=n_ixc_0+1
          ind_ixc(n_ixc_0)=4
          n_chi_0=n_chi_0+2
          ind_chi(3)=n_chi_0-1
          ind_chi(4)=n_chi_0
        endif
      endif
c ------- Exchange-Correlation Kernel ------------------------------
      if(chi_basis=='PB') call theta_xc(0)
      if(chi_basis=='PW') call theta_xc_PW(0)
      if(cht) then
c --------- Response_0 on TAU ---------------------------------------
        allocate(resp_lda_tau_mm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                           ndim3_k(me_k+1),n_chi_0))
        allocate(resp_lda_tau_mi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                           ndim3_tn,ndim3_k(me_k+1),n_chi_0))
        allocate(resp_lda_tau_ii(nn2,ndim3_tn,ndim3_k(me_k+1),n_chi_0))
        resp_lda_tau_mm=0.d0
        resp_lda_tau_mi=0.d0
        resp_lda_tau_ii=0.d0
        if(nspin==1) then
          call pol_0_spin_spin(1,1,resp_lda_tau_mm,resp_lda_tau_mi,
     &                         resp_lda_tau_ii)
          resp_lda_tau_mm=resp_lda_tau_mm*4.d0/nrel
          resp_lda_tau_mi=resp_lda_tau_mi*4.d0/nrel
          resp_lda_tau_ii=resp_lda_tau_ii*4.d0/nrel
        else if(nspin==2) then
          if(ind_chi(1)/=0) then
c ------- A(upup+dndn) and A(upup-dndn) -----------------------------
            i=ind_chi(1)
            j=ind_chi(2)
            call pol_0_spin_spin(2,2,resp_lda_tau_mm(1,1,1,i),
     &                           resp_lda_tau_mi(1,1,1,1,i),
     &                           resp_lda_tau_ii(1,1,1,i)) ! up_up
            call pol_0_spin_spin(1,1,resp_lda_tau_mm(1,1,1,j),
     &                           resp_lda_tau_mi(1,1,1,1,j),
     &                           resp_lda_tau_ii(1,1,1,j)) ! dn_dn
c ------------- Form proper combinations --------------------------
            do ind_k=1,ndim3_k(me_k+1)
              do ind_tau=1,ndim3_tau
                do i1=1,ncmpl*ndim_pbmt(me_b+1)
                  tmp=resp_lda_tau_mm(i1,ind_tau,ind_k,i)
     &	            +resp_lda_tau_mm(i1,ind_tau,ind_k,j)
                  resp_lda_tau_mm(i1,ind_tau,ind_k,j)=
     &	            (resp_lda_tau_mm(i1,ind_tau,ind_k,i)
     &	            -resp_lda_tau_mm(i1,ind_tau,ind_k,j))*2.d0
                  resp_lda_tau_mm(i1,ind_tau,ind_k,i)=tmp*2.d0
                enddo
                do j1=1,nplw_gw
                  do i1=1,ncmpl*nd_b_pbmt(me_b+1)
                    tmp=resp_lda_tau_mi(i1,j1,ind_tau,ind_k,i)
     &	               +resp_lda_tau_mi(i1,j1,ind_tau,ind_k,j)
                    resp_lda_tau_mi(i1,j1,ind_tau,ind_k,j)=
     &	            (resp_lda_tau_mi(i1,j1,ind_tau,ind_k,i)
     &	            -resp_lda_tau_mi(i1,j1,ind_tau,ind_k,j))*2.d0
                    resp_lda_tau_mi(i1,j1,ind_tau,ind_k,i)=tmp*2.d0
                  enddo
                enddo
                do i1=1,nn2
                  tmp=resp_lda_tau_ii(i1,ind_tau,ind_k,i)
     &	             +resp_lda_tau_ii(i1,ind_tau,ind_k,j)
                  resp_lda_tau_ii(i1,ind_tau,ind_k,j)=
     &	            (resp_lda_tau_ii(i1,ind_tau,ind_k,i)
     &	            -resp_lda_tau_ii(i1,ind_tau,ind_k,j))*2.d0
                  resp_lda_tau_ii(i1,ind_tau,ind_k,i)=tmp*2.d0
                enddo
              enddo
            enddo
          endif
          if(ind_chi(3)/=0) then
c ------- A(dnup+updn) and B(dnup-updn) -----------------------------
            i=ind_chi(3)
            j=ind_chi(4)
            call pol_0_spin_spin(1,2,resp_lda_tau_mm(1,1,1,i),
     &                           resp_lda_tau_mi(1,1,1,1,i),
     &                           resp_lda_tau_ii(1,1,1,i)) ! dn_up
            call pol_0_spin_spin(2,1,resp_lda_tau_mm(1,1,1,j),
     &                           resp_lda_tau_mi(1,1,1,1,j),
     &                           resp_lda_tau_ii(1,1,1,j)) ! up_dn
c ------------- Form proper A and B combinations for (+-) component --
            do ind_k=1,ndim3_k(me_k+1)
              do ind_tau=1,ndim3_tau
                do i1=1,ncmpl*ndim_pbmt(me_b+1)
                  tmp=resp_lda_tau_mm(i1,ind_tau,ind_k,i)
     &	            +resp_lda_tau_mm(i1,ind_tau,ind_k,j)
                  resp_lda_tau_mm(i1,ind_tau,ind_k,j)=
     &	            (resp_lda_tau_mm(i1,ind_tau,ind_k,i)
     &	            -resp_lda_tau_mm(i1,ind_tau,ind_k,j))*4.d0
                  resp_lda_tau_mm(i1,ind_tau,ind_k,i)=tmp*4.d0
                enddo
                do j1=1,nplw_gw
                  do i1=1,ncmpl*nd_b_pbmt(me_b+1)
                    tmp=resp_lda_tau_mi(i1,j1,ind_tau,ind_k,i)
     &	               +resp_lda_tau_mi(i1,j1,ind_tau,ind_k,j)
                    resp_lda_tau_mi(i1,j1,ind_tau,ind_k,j)=
     &	            (resp_lda_tau_mi(i1,j1,ind_tau,ind_k,i)
     &	            -resp_lda_tau_mi(i1,j1,ind_tau,ind_k,j))*4.d0
                    resp_lda_tau_mi(i1,j1,ind_tau,ind_k,i)=tmp*4.d0
                  enddo
                enddo
                do i1=1,nn2
                  tmp=resp_lda_tau_ii(i1,ind_tau,ind_k,i)
     &	             +resp_lda_tau_ii(i1,ind_tau,ind_k,j)
                  resp_lda_tau_ii(i1,ind_tau,ind_k,j)=
     &	            (resp_lda_tau_ii(i1,ind_tau,ind_k,i)
     &	            -resp_lda_tau_ii(i1,ind_tau,ind_k,j))*4.d0
                  resp_lda_tau_ii(i1,ind_tau,ind_k,i)=tmp*4.d0
                enddo
              enddo
            enddo
          endif
        endif
      endif
      call timel('**** RESP_LDA_PREPARE finished *****')
      end
