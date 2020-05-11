      subroutine vertex_pol_uniform(key,dpw)
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
	  real*8, intent(out) :: dpw
      logical :: q_present
      integer :: k,ispin,iq,i,it_vrt,n,ib,i_omega,i_tau,kst,it,
     &           ind_omega,kf,j,nrr3
      real*8 :: de,green0,lam_max,conv_1
      complex*16 :: cc
      integer, allocatable :: indrr3(:,:,:)
      real*8, allocatable :: v4_r(:,:,:),w4_r(:,:,:,:),
     &                       w_tau(:,:,:,:),g_tau(:,:,:,:,:,:),
     &                       g_om(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),vmi(:,:,:),
     &                       wmi(:,:,:,:),vim(:,:,:),
     &                       wim(:,:,:,:),vii(:),
     &                       wii(:,:),conv_table(:)
      complex*16, allocatable :: k0_pw(:,:,:,:,:),ax(:,:,:,:),
     &                           k_pw(:,:,:,:,:),pol(:),
     &                           lambda_dyn(:,:,:,:,:),
     &                           lambda_stat(:,:,:,:),v2_q(:,:,:),
     &                           z_red(:,:,:,:),
     &                           a_red(:,:,:,:),
     &                           vrt_dyn_old(:,:,:,:,:),
     &                           vrt_stat_old(:,:,:,:),
     &                           g_omega(:,:,:,:,:),phf(:,:)
      dpw=0.d0
      q_present=.false.
      if(nspin==2) q_present=.true.
      if(nspin==1.and.key==0) q_present=.true.
      n=nfun_red
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
c ---------------------------------------------------------------------  
      call norma_wf_red
c ---------- Interaction matrices in Q-orthonormal product basis ----
      allocate(w_tau(n_pbtot_red,n_pbtot_red,ndim3_tau,nqdiv_c))
      call w_tau_prepare_red(w_tau)
      if(q_present) call set_wnu_vrt
c ---------------------------------------------------------------------      
      allocate(v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c))
      v2_q=(0.d0,0.d0)
      do k=1,npnt_c
        kf=k_a_from_c(k)
	    do kst=1,k_star_c(k)
	      iq=k_list_c(kst,k)
	      v2_q(:,:,iq)=v_red_q(:,:,kf)
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
     &                          2*nbndf_bnd**2*(n_omega+1)*npnt_c*nspin,
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
	  allocate(lambda_dyn(nbndf_bnd,nbndf_bnd,0:n_omega,nqdiv_c,
     &                    nspin))
	  allocate(vrt_stat_old(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
	  allocate(vrt_dyn_old(nbndf_bnd,nbndf_bnd,0:n_omega,nqdiv_c,
     &                     nspin))
c --------------------------------------------------------------------
      allocate(conv_table(0:iter_ladder_xi))
      conv_table=0.d0
      allocate(pol(nspin))
c ------------------------------------------------------------------
	  pol=(0.d0,0.d0)
      lam_max=1.d0/sqrt(amega)
      conv_table(0)=lam_max
c --------------------------------------------------------------------
	  lambda_stat=(0.d0,0.d0)
	  lambda_dyn=(0.d0,0.d0)
c -------- Get zero order generalized susceptibility -------------------
      call k0_bnd_uniform(key,k0_pw,g_omega,g_x,ex)
c ----------------------------------------------------------------------
      vrt_stat_old=(0.d0,0.d0)
      vrt_dyn_old=(0.d0,0.d0)
      k_pw=(0.d0,0.d0)
      do it_vrt=1,iter_ladder_xi
        k_pw=k0_pw+k_pw
        call vertex_loop_uniform(q_present,it_vrt,k_pw,fif0,z_red,a_red,
     &                           phf,ax,v4_r,vmi,vim,vii,w4_r,wmi,wim,
     &                           wii,v2_q,lambda_stat,lambda_dyn,nspin,
     &                           nrr3,indrr3)
c --------------------------------------------------------------------
        conv_1=0.d0
        do ispin=1,nspin
          do k=1,nqdiv_c
            do i_omega=0,n_omega
              do j=1,nbndf_bnd
                do i=1,nbndf_bnd
                  cc=lambda_dyn(i,j,i_omega,k,ispin)
     &              +lambda_stat(i,j,k,ispin)
     &              -vrt_dyn_old(i,j,i_omega,k,ispin)
     &              -vrt_stat_old(i,j,k,ispin)
                  conv_1=max(conv_1,abs(cc))
                enddo
              enddo
            enddo
          enddo
        enddo
        conv_table(it_vrt)=
     &              max(conv_1,conv_table(it_vrt))
        vrt_dyn_old=lambda_dyn
        vrt_stat_old=lambda_stat
        do ispin=1,nspin
          do k=1,nqdiv_c
            if(it_vrt<iter_ladder_xi) then
c -------- Get correction to the generalized susceptibility --------
              call k_bnd_uniform(ispin,k,k_pw(1,1,0,k,ispin),
     &                          g_omega(1,1,0,1,ispin),
     &                          lambda_dyn(1,1,0,k,ispin),
     &                          lambda_stat(1,1,k,ispin),
     &                          g_x(1,0,1,ispin),ex(1,1,ispin))
            else
c -------------- Correction to the Polarizability ----------------------
              call dp_bnd_uniform(k,g_omega(1,1,0,1,ispin),
     &                           lambda_dyn(1,1,0,k,ispin),
     &                           pol(ispin),lambda_stat(1,1,k,ispin),
     &                           g_x(1,0,1,ispin),ex(1,1,ispin),
     &                           ispin)
            endif
          enddo  !! over k
        enddo  !! over ispin
      enddo  !! over it_vrt
      if(nspin==1) then
        dpw=dpw+pol(1)*2.d0
      else if(nspin==2) then
        if(key==0) then
          dpw=dpw+pol(1)+pol(2)
        else if(key==1) then
          dpw=dpw-pol(1)+pol(2)
        endif
      endif
      deallocate(pol)
c ------------------------------------------------------------------
	  if(maswrk) then
	    write(iun,'(a34,3f12.5)')' Convergency of Ladder eqn. '
	    do it=0,iter_ladder_xi
	      write(iun,'(i3,8e9.2)')it,conv_table(it)
	    enddo
      endif
      deallocate(v4_r,w4_r,v2_q,vmi,wmi,k0_pw,k_pw,lambda_dyn,ax,
     &           z_red,lambda_stat,g_x,ex,a_red,
     &           vim,wim,vii,wii,conv_table,vrt_dyn_old,vrt_stat_old,
     &           g_omega,phf,nind_fif,ind_fif,fif0,indrr3)
      if(q_present) deallocate(g_rs_mm,g_rs_mi,g_rs_im,g_rs_ii,wnu_vrt)
      call timel('*** VERTEX_POL_UNIFORM finished ****')
      end
