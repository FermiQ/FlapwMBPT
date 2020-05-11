      subroutine resp_lda(ch,n_ch0,ind_ch0,n_ix0,ind_ix0)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ch
      integer, intent(in) :: n_ch0,ind_ch0(n_ch0),n_ix0,ind_ix0(n_ix0)
      character*2 :: ch0
      integer :: ind_k,i,n_ch
      real*8, allocatable :: resp_stat(:,:),mtp(:)
      complex*16, allocatable :: xi0_nu(:,:,:),xi_nu(:,:,:),
     &                           eps1(:,:),chi_mt(:,:),eps_tns(:,:,:),
     &                           chi0_mt(:,:),chi_sort(:,:),
     &                           chi0_sort(:,:),chi_neig(:,:,:)
      if(ch=='___00') n_ch=1
      if(ch=='___ZZ') n_ch=1
      if(ch=='00_ZZ') n_ch=3
      if(ch=='___XY') n_ch=1
      if(ch=='___00'.or.ch=='00_ZZ') then
        if(.not.allocated(p_head_rsp))
     &               allocate(p_head_rsp(3,3,ndim3_nu))
        if(.not.allocated(p_wing_rsp))
     &               allocate(p_wing_rsp(3,n_pbmt+nplwgw(1),ndim3_nu))
      endif
      allocate(xi0_nu(0:n_nu,n_line,n_ch))
      allocate(xi_nu(0:n_nu,n_line,n_ch))
      allocate(eps1(0:n_nu,n_line))
      eps1=(0.d0,0.d0)
      allocate(eps_tns(3,3,0:n_nu))
      eps_tns=(0.d0,0.d0)
      xi0_nu=(0.d0,0.d0)
      xi_nu=(0.d0,0.d0)
      allocate(resp_stat(npnt,n_ch))
      resp_stat=0.d0
      allocate(chi_mt(natom,ndim3_nu))
      chi_mt=(0.d0,0.d0)
      allocate(chi0_mt(natom,ndim3_nu))
      chi0_mt=(0.d0,0.d0)
      allocate(chi_sort(0:n_nu,nsort))
      allocate(chi0_sort(0:n_nu,nsort))
      allocate(chi_neig(max_sites,nsort,ndim3_nu))
      chi_neig=(0.d0,0.d0)
      allocate(mtp(n_pbmt))
      mtp=0.d0
      if(me_k==0) mtp=pw_pb(1:n_pbmt,1)*sqrt(amega)
      if(nproc_k/=1) call dgop(mtp,n_pbmt,'  +',comm_k)
c -------------------------------------------------------------------
      do ind_k=1,ndim3_k(me_k+1)
        call ixc_q_nu_pb_0(ch,ind_k,xi0_nu,xi_nu,n_ch0,ind_ch0,n_ix0,
     &                     ind_ix0,n_ch,resp_stat,eps1,eps_tns,chi_mt,
     &                     chi0_mt,chi_neig,mtp)
      enddo  !! over ind_k
c -------------------------------------------------------------------
      if(nproc_t/=1) then
        call dgop(xi0_nu,2*(n_nu+1)*n_line*n_ch,'  +',comm_t)
        call dgop(xi_nu,2*(n_nu+1)*n_line*n_ch,'  +',comm_t)
        call dgop(resp_stat,npnt*n_ch,'  +',comm_t)
        call dgop(eps1,2*(n_nu+1)*n_line,'  +',comm_t)
        call dgop(eps_tns,18*(n_nu+1),'  +',comm_t)
      endif
      if(nproc_k/=1) then
        call dgop(xi0_nu,2*(n_nu+1)*n_line*n_ch,'  +',comm_k)
        call dgop(xi_nu,2*(n_nu+1)*n_line*n_ch,'  +',comm_k)
        call dgop(resp_stat,npnt*n_ch,'  +',comm_k)
        call dgop(eps1,2*(n_nu+1)*n_line,'  +',comm_k)
        call dgop(eps_tns,18*(n_nu+1),'  +',comm_k)
        call dgop(chi_mt,2*natom*ndim3_nu,'  +',comm_k)
        call dgop(chi0_mt,2*natom*ndim3_nu,'  +',comm_k)
        call dgop(chi_neig,2*max_sites*nsort*ndim3_nu,'  +',comm_k)
      endif
      xi_nu=xi_nu*amega
      xi0_nu=xi0_nu*amega
      resp_stat=resp_stat*amega
c -------------------------------------------------------------------
      if(ch=='___00') then
        call output_e_tensor(eps_tns,'_E')
        call output_boson_q_nu(eps1,'E1__00',n_line,'a')
        call output_boson_q_nu(xi0_nu,'X0__00',n_line,'a')
        call output_boson_q_nu(xi_nu,'X___00',n_line,'a')
        call output_xi_q('Rsp_00_STAT',resp_stat(1,1))
      endif
      if(ch=='___ZZ') then
        call output_boson_q_nu(xi0_nu,'X0__ZZ',n_line,'a')
        call output_boson_q_nu(xi_nu,'X___ZZ',n_line,'a')
        call output_xi_q('Rsp_ZZ_STAT',resp_stat(1,1))
        call output_chi_mt('Rsp_ZZ_MT__',chi_mt,chi0_mt,chi_sort,
     &                     chi0_sort)
        call output_chi_neig('Rsp_ZZ_MT__',chi_neig)
        call from_nu_to_rax_mt('ZZ',chi_sort)
        call from_nu_to_rax_mt('Z0',chi0_sort)
      endif
      if(ch=='00_ZZ') then
        call output_e_tensor(eps_tns,'_E')
        call output_boson_q_nu(eps1,'E1__00',n_line,'a')
        call output_boson_q_nu(xi0_nu(0,1,1),'X0__00',n_line,'a')
        call output_boson_q_nu(xi_nu(0,1,1),'X___00',n_line,'a')
        call output_boson_q_nu(xi0_nu(0,1,2),'X0__ZZ',n_line,'a')
        call output_boson_q_nu(xi_nu(0,1,2),'X___ZZ',n_line,'a')
        call output_boson_q_nu(xi0_nu(0,1,3),'X0__0Z',n_line,'a')
        call output_boson_q_nu(xi_nu(0,1,3),'X___0Z',n_line,'a')
        call output_xi_q('Rsp_00_STAT',resp_stat(1,1))
        call output_xi_q('Rsp_ZZ_STAT',resp_stat(1,2))
        call output_chi_mt('Rsp_ZZ_MT__',chi_mt,chi0_mt,chi_sort,
     &                     chi0_sort)
        call output_chi_neig('Rsp_ZZ_MT__',chi_neig)
        call from_nu_to_rax_mt('ZZ',chi_sort)
        call from_nu_to_rax_mt('Z0',chi0_sort)
      endif
      if(ch=='___XY') then
        call output_boson_q_nu(xi0_nu,'X0__+-',n_line,'a')
        call output_boson_q_nu(xi_nu,'X___+-',n_line,'a')
        call output_xi_q('Rsp_+-_STAT',resp_stat(1,1))
        call output_chi_mt('Rsp_+-_MT__',chi_mt,chi0_mt,chi_sort,
     &                     chi0_sort)
        call from_nu_to_rax_mt('+-',chi_sort)
        call from_nu_to_rax_mt('+0',chi0_sort)
      endif
c -------------- Real Axis --------------------------------------
      do i=1,n_ch
        if(ch=='___XY') then
	      call boson_real_axis('X0__+-',xi0_nu,n_line,.false.,1)
          xi0_nu=conjg(xi0_nu)
	      call boson_real_axis('X0__-+',xi0_nu,n_line,.false.,1)
          xi0_nu=conjg(xi0_nu)
	      call boson_real_axis('X___+-',xi_nu,n_line,.false.,1)
          xi_nu=conjg(xi_nu)
	      call boson_real_axis('X___-+',xi_nu,n_line,.false.,1)
          xi_nu=conjg(xi_nu)
        else
          if(ch=='___00') ch0='00'
          if(ch=='___ZZ') ch0='ZZ'
          if(ch=='00_ZZ'.and.i==1) ch0='00'
          if(ch=='00_ZZ'.and.i==2) ch0='ZZ'
          if(ch=='00_ZZ'.and.i==3) ch0='0Z'
	      call boson_real_axis('X0__'//ch0,xi0_nu(0,1,i),n_line,
     &                         .false.,2)
          call boson_real_axis('X___'//ch0,xi_nu(0,1,i),n_line,.false.,
     &                         2)
          if(ch=='___00'.or.ch=='00_ZZ') then
	        call boson_real_axis('E1____',eps1,n_line,.true.,2)
          endif
        endif
      enddo
      deallocate(resp_stat,xi0_nu,xi_nu,eps1,chi_mt,chi0_mt,eps_tns,
     &           chi_sort,chi0_sort,chi_neig)
      call timel('**** RESP_LDA finished *************')
      end
