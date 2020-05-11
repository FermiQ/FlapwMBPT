      subroutine vertex_pol_pb_static_gw(key,dii,dij)
c     key = 0   for 00-vertex
c     key = 1   for ZZ-vertex
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      complex*16, intent(out) :: dii(n_pbtot_red,n_pbtot_red,nc_line),
     &                           dij(n_pbtot_red,n_pbtot_red,nc_line)
      logical :: q_present
      integer :: k,ispin,iq,i,it_vrt,s2,np,n,ib,j,ind_omega,mx2,i_tau,
     &           kq,kst,iq_line,it,i_omega,ind_s2,kf,ind_qcl,iqf,ka,
     &           nrr3
      real*8 :: de,green0,v(3),gtild(3,2),lam_max,conv_1
      complex*16 :: cc
      integer, allocatable :: nd_np(:),n_np(:),ndim_kcl(:),n_mpi_kcl(:),
     &                        indrr3(:,:,:)
      real*8, allocatable :: v4_r(:,:,:),w4_r(:,:,:,:),
     &                       g_tau(:,:,:,:,:,:),w_tau(:,:,:,:),
     &                       g_om(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),vmi(:,:,:),
     &                       wmi(:,:,:,:),vim(:,:,:),
     &                       wim(:,:,:,:),vii(:),fif_full(:,:,:,:),
     &                       wii(:,:),conv_table(:,:)
      complex*16, allocatable :: k0_pw(:,:,:,:,:),
     &                           k_pw(:,:,:,:,:),pol(:,:,:),
     &                           lambda_dyn(:,:,:,:,:,:),
     &                           lambda_stat(:,:,:,:),v2_q(:,:,:),
     &                           z_red(:,:,:,:),ps3(:,:,:,:,:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           vrt_dyn_old(:,:,:,:,:,:),
     &                           vrt_stat_old(:,:,:,:),phf(:,:),
     &                           g_omega(:,:,:,:,:),
     &                           z_full(:,:,:,:),a_full(:,:,:,:)
      dii=(0.d0,0.d0)
      dij=(0.d0,0.d0)
      q_present=.false.
      if(nspin==2) q_present=.true.
      if(nspin==1.and.key==0) q_present=.true.
      n=nfun_red
      mx2=maxel_red**2
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
c ---------------------------------------------------------------------  
      call norma_wf_red
c ------ Here for zero-order vertex in full basis ------------------      
      allocate(fif_full(maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ---------- Interaction matrices in Q-orthonormal product basis ----
      allocate(w_tau(n_pbtot_red,n_pbtot_red,ndim3_tau,nqdiv_c))
      call w_tau_prepare_red(w_tau)
      if(q_present) call set_wnu_vrt
c ---------------------------------------------------------------------    
      allocate(v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c))
      v2_q=(0.d0,0.d0)
      do k=1,npnt_c
        ka=k_a_from_c(k)
	    do kst=1,k_star_c(k)
	      iq=k_list_c(kst,k)
	      v2_q(:,:,iq)=v_red_q(:,:,ka)
	      call sym_w_red(iq,k,v2_q(1,1,iq))
        enddo   !! over kst
      enddo   !! over k
c --------- Interaction matrices in 2-site (MM) form ------------------
      allocate(v4_r(n_pbmt_red,n_pbmt_red,nqdiv_c))
      allocate(w4_r(n_pbmt_red,n_pbmt_red,nqdiv_c,0:n_tau/2))
      call w4_for_bnd_mm(v4_r,w4_r,w_tau,v2_q,n_pbtot_red)
c --------- Interaction matrices in (MI) form ------------------
      allocate(vmi(n_pbmt_red,nr_full_red,nqdiv_c))
      allocate(wmi(n_pbmt_red,nr_full_red,nqdiv_c,0:n_tau/2))
      call w4_for_bnd_mi(vmi,wmi,w_tau,v2_q,phf)
c --------- Interaction matrices in (IM) form ------------------
      allocate(vim(nr_full_red,n_pbmt_red,nqdiv_c))
      allocate(wim(nr_full_red,n_pbmt_red,nqdiv_c,0:n_tau/2))
      call w4_for_bnd_im(vim,wim,w_tau,v2_q,phf)
c --------- Interaction matrices in (II) form ------------------
      nrr3=nr_red_red*nr_full_red*nqdiv_c
      allocate(indrr3(nr_full_red,nr_full_red,nqdiv_c))
      allocate(vii(nrr3))
      allocate(wii(nrr3,0:n_tau/2))
      call w4_for_bnd_ii(vii,wii,w_tau,v2_q,phf,nrr3,indrr3)
      deallocate(w_tau)
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_om(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_om)
c -----------------
	  allocate(g_omega(nbndf_bnd,nbndf_bnd,0:n_omega,npnt_c,nspin))
	  g_omega=(0.d0,0.d0)
	  do ispin=1,nspin
	    do k=1,npnt_c
          kf=k_a_from_c(k)
	      n=n_low_bnd(kf,ispin)
	      do ind_omega=1,ndim3_omega
            i_omega=me_t*ndim3_omega+ind_omega-1
	        call ferm_unpack_omega(g_omega(1,1,i_omega,k,ispin),
     &	                           g_om(1,1,1,ind_omega,k,ispin),n,
     &                             nbndf_bnd,nbndf_bnd)
          enddo
        enddo
      enddo
      deallocate(g_om)
	  if(nproc_t/=1) call dgop(g_omega,
     &                        2*nbndf_bnd**2*(n_omega+1)*npnt_c*nspin,
     &                           '  +',comm_t)
c -----------------
      if(q_present) then
        allocate(g_rs_mm(nrel,nfun_red,nfun_red,0:n_tau,nqdiv_c,nspin))
        allocate(g_rs_mi(nrel,nfun_red,nrel*nr_full_red,0:n_tau,nqdiv_c,
     &                   nspin))
        allocate(g_rs_im(nrel,nrel*nr_full_red,nfun_red,0:n_tau,nqdiv_c,
     &                   nspin))
        allocate(g_rs_ii(nrel,nrel*nr_red_red,nrel*nr_full_red,0:n_tau,
     &                   nqdiv_c,nspin))
        do ispin=1,nspin
          call g_real_space(g_tau(1,1,1,1,1,ispin),
     &                      g_rs_mm(1,1,1,0,1,ispin),
     &                      g_rs_mi(1,1,1,0,1,ispin),
     &                      g_rs_im(1,1,1,0,1,ispin),
     &                      g_rs_ii(1,1,1,0,1,ispin),z_red(1,1,1,ispin),
     &                      a_red(1,1,1,ispin),phf,ispin)
        enddo
        call timel('**** G_REAL_SPACE finished *********')
      endif
      deallocate(g_tau)
c ---- Exchange eigen values relative to Chem_POT ------------------
      allocate(ex(nbndf_bnd,npnt_c,nspin))
	  do ispin=1,nspin
	    do k=1,npnt_c
          kf=k_a_from_c(k)
     	  do i=1,n_low_bnd(kf,ispin)
     	    ib=ind_bands_bnd(i,kf,ispin)
            ex(i,k,ispin)=e_bnd(ib,kf,ispin)-chem_pot
          enddo
	    enddo   !! over k
	  enddo  !! over ispin
c --------- Exchange Green's function ------------------------------
      allocate(g_x(nbndf_bnd,0:n_tau,npnt_c,nspin))
	  do ispin=1,nspin
	    do k=1,npnt_c
          kf=k_a_from_c(k)
     	  do i=1,n_low_bnd(kf,ispin)
            de=ex(i,k,ispin)
            do i_tau=0,n_tau
	          g_x(i,i_tau,k,ispin)=green0(de,tau_mesh(i_tau))
            enddo
          enddo
	    enddo   !! over k
	  enddo  !! over ispin
c --------------------------------------------------------------------
	  allocate(k0_pw(nbndf_bnd,nbndf_bnd,0:n_tau,nqdiv_c,nspin))
	  allocate(k_pw(nbndf_bnd,nbndf_bnd,0:n_tau,nqdiv_c,nspin))
	  allocate(lambda_stat(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
	  allocate(lambda_dyn(nbndf_bnd,nbndf_bnd,0:n_omega,2,nqdiv_c,
     &                    nspin))
	  allocate(vrt_stat_old(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
	  allocate(vrt_dyn_old(nbndf_bnd,nbndf_bnd,0:n_omega,2,nqdiv_c,
     &                     nspin))
c --------------------------------------------------------------------
      allocate(conv_table(0:iter_ladder_xi,nc_line))
      conv_table=0.d0
      allocate(pol(n_pbtot_red,n_pbtot_red,nspin))
      pol=(0.d0,0.d0)
      allocate(ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,nqdiv_c,nspin))
c ------------------------------------------------------------------ 
      allocate(nd_np(nproc_t))
      allocate(n_np(nproc_t))
	  allocate(ndim_kcl(nproc_k))
	  allocate(n_mpi_kcl(nproc_k))
	  call size_shift_par(nc_line,nproc_k,ndim_kcl,n_mpi_kcl)
	  do ind_qcl=1,ndim_kcl(me_k+1)
        iq_line=n_mpi_kcl(me_k+1)+ind_qcl
	    iq=kline_in_npnt_c(iq_line)
        iqf=k_a_from_c(iq)
        ka=k_a_from_c(iq)
        np=n_pbmt_red+nplwgw_red(ka)
        call size_shift_par(np,nproc_t,nd_np,n_np)
        ps3=(0.d0,0.d0)
        do ispin=1,nspin
          do k=1,nqdiv_c
            v=pnt_c(:,k)-pnt_c(:,iq)
            call zone1_number(v,rb0_c,ndiv_c,kq)
            kq=index_k1_c(kq)    !! for K-Q
            gtild(:,1)=0.d0
            gtild(:,2)=v-pnt_c(:,kq)
            call ppm_factors_full(k,kq,iq,fif_full(1,1,1,ispin),
     &	                          z_full(1,1,1,ispin),
     &                            a_full(1,1,1,ispin),
     &                            ps3(1,1,1,k,ispin),gtild,
     &                            n_low_bnd(1,ispin))
          enddo   !! over k
        enddo   !! over ispin
        lam_max=maxval(abs(ps3))
        conv_table(0,iq_line)=lam_max
        pol=(0.d0,0.d0)
        do ind_s2=1,nd_np(me_t+1)
          s2=n_np(me_t+1)+ind_s2
c --------------------------------------------------------------------
          lambda_stat=(0.d0,0.d0)
          lambda_dyn=(0.d0,0.d0)
c -------- Get zero order generalized susceptibility -------------------
          call k0_bnd_static(key,iq,k0_pw,g_omega,ps3,s2,g_x,ex)
c ----------------------------------------------------------------------
          vrt_stat_old=(0.d0,0.d0)
          vrt_dyn_old=(0.d0,0.d0)
          k_pw=(0.d0,0.d0)
          do it_vrt=1,iter_ladder_xi
            k_pw=k0_pw+k_pw
            call vertex_bnd_loop(q_present,iq,1,it_vrt,1,s2,n_omega+1,1,
     &                           k_pw,fif0,z_red,a_red,phf,ax,
     &                           v4_r,vmi,vim,vii,w4_r,wmi,wim,wii,v2_q,
     &                           lambda_stat,lambda_dyn,nspin,0,nrr3,
     &                           indrr3)
c --------------------------------------------------------------------
            conv_1=0.d0
            do ispin=1,nspin
              do k=1,nqdiv_c
                do i_omega=0,n_omega
                  do it=1,2
                    do j=1,nbndf_bnd
                      do i=1,nbndf_bnd
                        cc=lambda_dyn(i,j,i_omega,it,k,ispin)
     &                    +lambda_stat(i,j,k,ispin)
     &                    -vrt_dyn_old(i,j,i_omega,it,k,ispin)
     &                    -vrt_stat_old(i,j,k,ispin)
                        conv_1=max(conv_1,abs(cc))
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
            conv_table(it_vrt,iq_line)=
     &              max(conv_1,conv_table(it_vrt,iq_line))
            vrt_dyn_old=lambda_dyn
            vrt_stat_old=lambda_stat
            do ispin=1,nspin
              do k=1,nqdiv_c
                if(it_vrt<iter_ladder_xi) then
c -------- Get correction to the generalized susceptibility --------
                  call k_bnd_static(ispin,k,iq,
     &                              k_pw(1,1,0,k,ispin),
     &                              g_omega(1,1,0,1,ispin),
     &                              lambda_dyn(1,1,0,1,k,ispin),
     &                              lambda_stat(1,1,k,ispin),
     &                              g_x(1,0,1,ispin),ex(1,1,ispin))
                else
c -------------- Correction to the Polarizability ----------------------
                  call dp_bnd_static(k,iq,g_omega(1,1,0,1,ispin),
     &                               lambda_dyn(1,1,0,1,k,ispin),
     &                               ps3(1,1,1,k,ispin),
     &                               pol(1,s2,ispin),
     &                               lambda_stat(1,1,k,ispin),
     &                               g_x(1,0,1,ispin),ex(1,1,ispin),
     &                               ispin,n_pbtot_red)
                endif
              enddo  !! over k
            enddo  !! over ispin
          enddo  !! over it_vrt
        enddo   !! over ind_s2
        if(nproc_t/=1) call dgop(pol,2*n_pbtot_red**2*nspin,'  +',
     &	                         comm_t)
        if(nspin==1) then
          dii(:,:,iq_line)=dii(:,:,iq_line)+pol(:,:,1)*2.d0
        else if(nspin==2) then
          if(key==0) then
            dii(:,:,iq_line)=dii(:,:,iq_line)+pol(:,:,1)+pol(:,:,2)
            dij(:,:,iq_line)=dij(:,:,iq_line)-pol(:,:,1)+pol(:,:,2)
          else if(key==1) then
            dii(:,:,iq_line)=dii(:,:,iq_line)-pol(:,:,1)+pol(:,:,2)
            dij(:,:,iq_line)=dij(:,:,iq_line)+pol(:,:,1)+pol(:,:,2)
          endif
        endif
      enddo  !! over ind_qcl
      deallocate(pol,ps3)
      if(nproc_k/=1) then
        call dgop(dii,2*n_pbtot_red**2*nc_line,'  +',comm_k)
        call dgop(dij,2*n_pbtot_red**2*nc_line,'  +',comm_k)
        call dgop(conv_table,nc_line*(iter_ladder_xi+1),'  +',comm_k)
      endif
      if(nproc_t/=1) call dgop(conv_table,nc_line*(iter_ladder_xi+1),
     &                           'max',comm_t)
c ------------------------------------------------------------------
	if(maswrk) then
	  do iq=1,nc_line
	    k=kline_in_npnt_c(iq)
	    write(iun,'(a34,3f12.5)')
     &	      ' Convergency of Ladder eqn. for Q= ',
     &	              (pnt_c(i,k),i=1,3)
	    write(iun,'(a5,8(i3,6x))')'It\Nu',0,1,2,3,4,5,6,7
	    do it=0,iter_ladder_xi
	      write(iun,'(i3,8e9.2)')it,conv_table(it,iq)
	    enddo
        enddo
      endif
c ------------------------------------------------------------------
      deallocate(v4_r,w4_r,v2_q,vmi,wmi,k0_pw,k_pw,lambda_dyn,ax,
     &           z_red,lambda_stat,g_x,ex,a_red,indrr3,
     &           vim,wim,vii,wii,conv_table,vrt_dyn_old,vrt_stat_old,
     &           g_omega,nd_np,phf,nind_fif,ind_fif,fif0,
     &           n_np,fif_full,z_full,a_full,ndim_kcl,n_mpi_kcl)
      if(q_present) deallocate(g_rs_mm,g_rs_mi,g_rs_im,g_rs_ii,wnu_vrt)
      call timel(' VERTEX_POL_PB_STATIC_GW finished **')
      end
