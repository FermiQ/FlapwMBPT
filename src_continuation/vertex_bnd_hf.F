      subroutine vertex_bnd_hf(iter_w,appr)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*3, intent(in) :: appr
      integer, intent(in) :: iter_w
      integer :: k,ispin,iq,ind_nu,i_nu,i,it_vrt,s2,np,nit,i0,iatom,
     &           ind_q,ib,mx2,i_tau,kq,kst,iq_line,ii,isort,khbr,
     &           ind_s2,idum,kf,ka,j,ind_k,n_pbt
      real*8 :: p_herm,check_hermit,ph,de,green0,v(3),gtild(3,2),
     &          lam_max,conv_1,pi2,vii0,phase
      complex*16 :: cc
      real*8, allocatable :: v4_r(:,:,:),poli(:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),vmi(:,:,:),vim(:,:,:),
     &                       vii(:,:,:),conv_table(:,:,:),
     &                       fif_full(:,:,:,:,:),vloc(:,:,:),
     &                       fxc_mt(:,:,:,:)
      complex*16, allocatable :: k0_pw(:,:,:,:),k_pw(:,:,:,:),
     &                           pol(:,:,:,:),lambda_stat(:,:,:,:),
     &                           v2_q(:,:,:),z_red(:,:,:,:),ax(:,:,:,:),
     &                           ps3(:,:,:,:,:),a_red(:,:,:,:),
     &                           vrt_stat_old(:,:,:,:),z_full(:,:,:,:),
     &                           a_full(:,:,:,:),phf(:,:),v_tmp(:,:)
      pi2=pi+pi
      mx2=maxel_red**2
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
c ------------------------------------------------------------------
      call norma_wf_red
c ------ Here for zero-order vertex in full basis ------------------
      allocate(fif_full(nrel,maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ------------------------------------------------------------------
      khbr=0
      if(hybrid_type>2) khbr=1
      if(khbr==1) then
        allocate(v_red_b(n_pbtot_red,n_pbtot_red,npnt_c))
        v_red_b=(0.d0,0.d0)
        do ind_k=1,ndim_k_red(me_k+1)
          k=n_mpi_k_red(me_k+1)+ind_k
          ka=k_a_from_c(k)
          n_pbt=n_pbmt+nplwgw(ka)
          allocate(v_tmp(n_pbt,n_pbt))
c          call v_full_WWWWWWWWWW(v_tmp,nplwgw(ka),k,1,2,1)
          call rpb_from_fpb_full(v_red_b(1,1,k),v_tmp,ka,
     &                           n_pbtot_red,n_pbt)
          deallocate(v_tmp)
        enddo
        if(nproc_k/=1) call dgop(v_red_b,2*n_pbtot_red**2*npnt_c,'  +',
     &                           comm_k)
      endif
c ------------------------------------------------------------------
      if(appr=='wmt') then
        allocate(vloc(n_pbmtm_red,n_pbmtm_red,natom))
        call v_loc_from_red(vloc,khbr)
        if(key_hybrid==1) then
          vloc=alpha_hbr*vloc
          n_ixc_0=1
          ind_ixc(1)=1
          allocate(theta_mt(maxmt/nspin,n_ixc_0))
          call theta_mt_xc(0)
          allocate(fxc_mt(n_pbmtm_red,n_pbmtm_red,natom,n_ixc_0))
          call theta_mtm_red(fxc_mt)
          deallocate(theta_mt)
          vloc=vloc-fxc_mt(:,:,:,1)
          deallocate(fxc_mt)
        endif
      endif
c ------------------------------------------------------------------
      if(appr=='bnd') then
        allocate(v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c))
        v2_q=(0.d0,0.d0)
        do k=1,npnt_c
          ka=k_a_from_c(k)
          do kst=1,k_star_c(k)
            iq=k_list_c(kst,k)
            if(khbr==0) v2_q(:,:,iq)=v_red_q(:,:,ka)
            if(khbr==1) v2_q(:,:,iq)=v_red_b(:,:,k)
            call sym_w_red(iq,k,v2_q(1,1,iq))
          enddo   !! over kst
        enddo   !! over k
        if(key_hybrid==1) v2_q=alpha_hbr*v2_q
        call timel('**** VERTEX_BND_HF : v2_q finished *')
c --------- Interaction matrices in 2-site (MM) form ---------------
        allocate(v4_r(n_pbmt_red,n_pbmt_red,nqdiv_c))
        call v4_for_bnd_mm(v4_r,v2_q,n_pbtot_red)
c --------- Interaction matrices in (MI) form ------------------
        allocate(vmi(n_pbmt_red,nr_full_red,nqdiv_c))
        call v4_for_bnd_mi(vmi,v2_q,phf)
c --------- Interaction matrices in (IM) form ------------------
        allocate(vim(nr_full_red,n_pbmt_red,nqdiv_c))
        call v4_for_bnd_im(vim,v2_q,phf)
c --------- Interaction matrices in (II) form ------------------
        allocate(vii(nr_red_red,nr_full_red,nqdiv_c))
        call v4_for_bnd_ii(vii,v2_q,phf)
        deallocate(v2_q)
        if(key_hybrid==1) then
          n_ixc_0=1
          ind_ixc(1)=1
          allocate(theta_mt(maxmt/nspin,n_ixc_0))
          call theta_mt_xc(0)
          allocate(fxc_mt(n_pbmtm_red,n_pbmtm_red,natom,n_ixc_0))
          call theta_mtm_red(fxc_mt)
          deallocate(theta_mt)
          do iatom=1,natom
            isort=is(iatom)
            i0=iopb_red(iatom)
            do j=1,n_pbmt0_red(isort)
              do i=1,n_pbmt0_red(isort)
                v4_r(i0+i-1,i0+j-1,1)=v4_r(i0+i-1,i0+j-1,1)
     &                               -fxc_mt(i,j,iatom,1)
              enddo
            enddo
          enddo
          deallocate(fxc_mt)
          allocate(theta_int(nplwro,n_ixc_0))
          call theta_int_xc(0)
          do i0=1,nr_red_red
            i=r_pnt_ind_red(i0)
            vii0=0.d0
            do ii=1,nplwro
              phase=pi2*dot_product(gbs(:,ii),r_pnt_red(:,i))
              vii0=vii0+theta_int(ii,1)*dcmplx(cos(phase),sin(phase))
            enddo
            vii(i0,i,1)=vii(i0,i,1)-vii0
          enddo
        endif
      endif
      if(khbr==1) deallocate(v_red_b)
      call timel('*** VERTEX_BND_HF : all V finished *')
c ---- Exchange eigen values relative to Chem_POT ------------------
      allocate(ex(nbndf_bnd,npnt_c,nspin))
      ex=0.d0
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
      g_x=0.d0
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
c -------------
      allocate(k0_pw(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      allocate(k_pw(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      allocate(lambda_stat(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      allocate(vrt_stat_old(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
c --------------------------------------------------------------------
      allocate(conv_table(0:iter_w,0:7,nc_line))
      conv_table=0.d0
      allocate(pol(n_pbtot_red,n_pbtot_red,ndim3_nu,
     &             ndim_k_red(me_k+1)))
      pol=(0.d0,0.d0)
      allocate(ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,nqdiv_c,nspin))
      call timel('*** most allocations finished ******')
c ------------------------------------------------------------------
      call timel('**** VERTEX_BND_HF : loop starts ***')
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        ka=k_a_from_c(iq)
        iq_line=k_npnt_in_line_c(iq)
        np=n_pbmt_red+nplwgw_red(ka)
        ps3=(0.d0,0.d0)
        do ispin=1,nspin
          do k=1,nqdiv_c
            v=pnt_c(:,k)-pnt_c(:,iq)
            call zone1_number(v,rb0_c,ndiv_c,kq)
            kq=index_k1_c(kq)    !! for K-Q
            gtild(:,1)=0.d0
            gtild(:,2)=v-pnt_c(:,kq)
            call ppm_factors_full(k,kq,iq,fif_full(1,1,1,1,ispin),
     &                            z_full(1,1,1,ispin),
     &                            a_full(1,1,1,ispin),
     &                            ps3(1,1,1,k,ispin),gtild,
     &                            n_low_bnd(1,ispin))
          enddo   !! over k
        enddo   !! over ispin
        lam_max=maxval(abs(ps3))
        if(iq_line/=0) conv_table(0,:,iq_line)=lam_max/nproc_t
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          nit=iter_w
c --------------------------------------------------------------------
          do ind_s2=1,ndim4_pbr(me_b+1,iq)
            s2=n4_mpi_pbr(me_b+1,iq)+ind_s2
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND starts *****************')
c -------- Get zero order generalized susceptibility ---------------
            call k0_bnd_hf(i_nu,iq,k0_pw,ps3,s2,g_x,ex)
            if(iq==1.and.i_nu==0.and.s2==1)
     &       call timel('**** K0_BND finished ***************')
c ----------------------------------------------------------------------
            k_pw=(0.d0,0.d0)
            vrt_stat_old=(0.d0,0.d0)
            do it_vrt=1,nit
              if(iq==1.and.i_nu==0.and.s2==1.and.it_vrt==1)
     &          call timel('****** Iter-n starts ***************')
              k_pw=k0_pw+k_pw
              if(appr=='bnd') then
                lambda_stat=(0.d0,0.d0)
                call vertex_vk(iq,k_pw,z_red,a_red,phf,ax,v4_r,vmi,vim,
     &                         vii,lambda_stat,nspin,0,fif0,khbr)
              else if(appr=='wmt') then
                lambda_stat=(0.d0,0.d0)
                call vertex_vk_wmt(iq,k_pw,z_red,vloc,lambda_stat,nspin,
     &                             fif0)
              endif
c --------------------------------------------------------------------
              conv_1=0.d0
              do ispin=1,nspin
                do k=1,nqdiv_c
                  do j=1,nbndf_bnd
                    do i=1,nbndf_bnd
                      cc=lambda_stat(i,j,k,ispin)
     &                  -vrt_stat_old(i,j,k,ispin)
                      conv_1=max(conv_1,abs(cc))
                    enddo
                  enddo
                enddo
              enddo
              if(i_nu<=7.and.iq_line/=0) then
                conv_table(it_vrt,i_nu,iq_line)=
     &              max(conv_1,conv_table(it_vrt,i_nu,iq_line))
              endif
              vrt_stat_old=lambda_stat
              do ispin=1,nspin
                if(it_vrt<nit) then
                  k_pw(:,:,:,ispin)=(0.d0,0.d0)
                  do k=1,nqdiv_c
c -------- Get correction to the generalized susceptibility --------
                    call k_bnd_hf(ispin,i_nu,k,iq,k_pw(1,1,k,ispin),
     &                            lambda_stat(1,1,k,ispin),
     &                            g_x(1,0,1,ispin),ex(1,1,ispin))
                  enddo  !! over k
                else if(it_vrt==nit) then
c -------------- Correction to the Polarizability ------------------
                  do k=1,nqdiv_c
                    call dp_bnd_hf(i_nu,k,iq,lambda_stat(1,1,k,ispin),
     &                             ps3(1,1,1,k,ispin),
     &                             pol(1,s2,ind_nu,ind_q),ispin,
     &                             g_x(1,0,1,ispin),ex(1,1,ispin),
     &                             n_pbtot_red)
                  enddo  !! over k
                endif
              enddo  !! over ispin
            enddo  !! over it_vrt
          enddo   !! over ind_s2
          if(iq==1.and.i_nu==0)
     &       call timel('**** S2-loop finished **************')
        enddo   !! over ind_nu
        if(iq==1) call timel('**** Nu-loop finished **************')
        if(iq==1) call timel('**** END of Q=1 loop ***************')
      enddo  !! over ind_q
      call timel('**** END of Q loop *****************')
      deallocate(k0_pw,k_pw,lambda_stat,g_x,ex,vrt_stat_old,fif_full,
     &           z_full,a_full)
      if(appr=='bnd') deallocate(v4_r,vmi,vim,vii)
      if(appr=='wmt') deallocate(vloc)
c ---------------------------------------------------------------
      if(nproc_b/=1) call dgop(pol,
     &                2*n_pbtot_red**2*ndim3_nu*ndim_k_red(me_k+1),
     &                         '  +',comm_b)
      if(nspin==1.and.irel/=2) pol=2.d0*pol
      p_herm=0.d0
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        ka=k_a_from_c(iq)
        np=n_pbmt_red+nplwgw_red(ka)
        allocate(poli(np,np))
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          ph=check_hermit(pol(1,1,ind_nu,ind_q),np,n_pbtot_red)
          p_herm=max(p_herm,ph)
          call pack_hermit(pol(1,1,ind_nu,ind_q),poli,np,
     &                     n_pbtot_red,np,0.d0,1.d0)
          p_vertex_nu_q(1:np,1:np,ind_nu,ind_q)=poli
        enddo
        deallocate(poli)
      enddo
      if(nproc_t/=1) call dgop(p_herm,1,'max',comm_t)
      if(nproc_k/=1) call dgop(p_herm,1,'max',comm_k)
      if(maswrk) write(iun,*)' VERTEX_BND_HF : P_HERM = ',p_herm
      call output_p_bnd_nu(2,pol,n_pbtot_red)
      deallocate(ps3,z_red,a_red,phf,nind_fif,ind_fif,fif0,ax)
c ------------------------------------------------------------------
      idum=(iter_w+1)*8*nc_line
      if(nproc_b/=1) call dgop(conv_table,idum,'max',comm_b)
      if(nproc_k/=1) call dgop(conv_table,idum,'  +',comm_k)
      if(nproc_t/=1) call dgop(conv_table,idum,'  +',comm_t)
      call output_bs_conv(conv_table,iter_w)
      deallocate(pol,conv_table)
      call timel('**** VERTEX_BND_HF finished ********')
      end
