      subroutine vertex_p_loc(iter_w,iter_t)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iter_w,iter_t
      logical :: spin_flips
      integer :: ispin,ind_nu,i_nu,i,nom,it_vrt,s2,i_omega,np,nit,
     &           ind_q,nbb,ind_omega,mx2,i_tau,it,
     &           nomnu,idum,j,iatom,ind,ind_tau,irl,
     &           isort,n
      real*8 :: p_herm,check_hermit,ph,lam_max,conv_1,pi2
      complex*16 :: cc
      real*8, allocatable :: poli(:,:),gexa(:,:,:,:,:),g0(:,:),
     &                       gspl_geom(:,:,:,:,:,:),
     &                       gspl_asy(:,:,:,:,:,:),
     &                       gasy(:,:,:,:),fif(:,:,:,:,:),
     &                       w_tau(:,:,:),g_tau(:,:,:,:,:,:,:),
     &                       fifim(:,:,:,:,:),conv_table(:,:,:),
     &                       wloc(:,:,:,:),gtau(:,:,:,:,:),
     &                       vloc(:,:,:),wt(:,:,:)
      complex*16, allocatable :: k0_pw(:,:,:,:,:),trn(:,:),
     &                           k_pw(:,:,:,:,:),pol(:,:,:),
     &                           lambda_dyn(:,:,:,:,:),
     &                           lambda_stat(:,:,:),
     &                           g_om_nu(:,:,:,:,:),
     &                           z_red(:,:,:,:),
     &                           vrt_dyn_old(:,:,:,:,:),
     &                           vrt_stat_old(:,:,:),
     &                           g_omega(:,:,:,:,:)
      pi2=pi+pi
      mx2=maxel_red**2
      allocate(fif(nrel,mx2,n_pbmtm_red,nsort,nspin))
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_loc(fif,z_red)
      deallocate(fif)
      allocate(wloc(n_pbmtm_red,n_pbmtm_red,0:n_nu,natom))
      allocate(vloc(n_pbmtm_red,n_pbmtm_red,natom))
      call w_loc_from_red(vloc,wloc)
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nrel,maxel_red,maxel_red,2,ndim3_tau,nspin,natom))
      allocate(g_omega(maxel_red,maxel_red,ndim3_omega,nspin,natom))
      call g_loc_from_full(g_omega,g_tau,z_red)
c --------------------------------------------------------------------
      allocate(conv_table(0:iter_w,0:7,natom))
      conv_table=0.d0
c ------------------------------------------------------------------
      allocate(trn(n_omega1_max,0:n_nu))
      trn=(0.d0,0.d0)
      if(me_t==0) then
        do i_nu=0,n_nu
          trn(:,i_nu)=tau_from_omega_nu(:,1,1,i_nu)
        enddo
      endif
      if(nproc_t/=1) call dgop(trn,2*n_omega1_max*(n_nu+1),'  +',
     &                           comm_t)
c ---------------------------------------------------------------
      allocate(fifim(nrel,maxel_red,maxel_red,n_pbmtm_red,nsort))
      do iatom=1,natom
        isort=is(iatom)
        n=lfunm_red(isort)
        if(n==0) cycle
        np=n_pbmt0_red(isort)
        allocate(fif(nrel,n,n,np,nspin))
c ------------------------------------------------------------------
        do ispin=1,nspin
          call fi0_full_red(fifim,ispin,ispin)
          do j=1,n
            do i=1,n
              do irl=1,nrel
                fif(irl,i,j,1:np,ispin)=fifim(irl,i,j,1:np,isort)
              enddo
            enddo
          enddo
        enddo
        allocate(w_tau(np,np,0:n_tau/2))
        allocate(wt(n_pbmtm_red,n_pbmtm_red,0:n_tau/2))
        call nu_to_tau_spl(wloc(1,1,0,iatom),wt,n_pbmtm_red**2)
        do i_tau=0,n_tau/2
          w_tau(:,:,i_tau)=wt(1:np,1:np,i_tau)/2
        enddo
        deallocate(wt)
c -------------
        allocate(gtau(nrel,n,n,0:n_tau,nspin))
        gtau=0.d0
        do ispin=1,nspin
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=me_t*ndim3_tau+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
              gtau(:,:,:,i_tau,ispin)=
     &          g_tau(:,1:n,1:n,it,ind_tau,ispin,iatom)
            enddo
          enddo
        enddo
        if(nproc_t/=1) call dgop(gtau,n*n*nrel*(n_tau+1)*nspin,'  +',
     &                             comm_t)
        allocate(gexa(n,n,2,0:n_omega_exa,nspin))
        allocate(gspl_geom(n,n,2,0:n_omega_geom+2,4,nspin))
        allocate(gspl_asy(n,n,2,n_omega_asy+1,4,nspin))
        allocate(gasy(n,n,2,nspin))
        gexa=0.d0
        gspl_geom=0.d0
        gspl_asy=0.d0
        gasy=0.d0
c --------- For Interpolation of G(w') ---------------------------------
        nbb=n**2
        allocate(g0(0:n_omega,2))
        do ispin=1,nspin
          do j=1,n
            do i=1,n
              g0=0.d0
              do ind_omega=1,ndim3_omega
                i_omega=me_t*ndim3_omega+ind_omega-1
                g0(i_omega,1)=real(g_omega(i,j,ind_omega,ispin,iatom))
                g0(i_omega,2)=imag(g_omega(i,j,ind_omega,ispin,iatom))
              enddo
              if(nproc_t/=1) call dgop(g0,2*(n_omega+1),'  +',
     &                                   comm_t)
              do it=1,2
                call spline_inhmg(w_omega(n_omega_exa-1),
     &                            g0(n_omega_exa-1:,it),
     &                            gspl_geom(i,j,it,:,1,ispin),
     &                            gspl_geom(i,j,it,:,2,ispin),
     &                            gspl_geom(i,j,it,:,3,ispin),
     &                            gspl_geom(i,j,it,:,4,ispin),
     &                            n_omega_geom+3,0,0.d0,0.d0)
                call spline_inhmg(xm_omega(0),
     &                            g0(n_omega_exa+n_omega_geom:,it),
     &                            gspl_asy(i,j,it,:,1,ispin),
     &                            gspl_asy(i,j,it,:,2,ispin),
     &                            gspl_asy(i,j,it,:,3,ispin),
     &                            gspl_asy(i,j,it,:,4,ispin),
     &                            n_omega_asy+1,0,0.d0,0.d0)
                gasy(i,j,it,ispin)=g0(n_omega,it)
                gexa(i,j,it,:,ispin)=g0(0:n_omega_exa,it)
              enddo
            enddo
          enddo
        enddo
        deallocate(g0)
        allocate(k0_pw(n,n,0:n_tau,2,nspin))
        allocate(k_pw(n,n,0:n_tau,2,nspin))
        allocate(lambda_stat(n,n,nspin))
        allocate(vrt_stat_old(n,n,nspin))
        allocate(pol(np,np,ndim3_nu))
        pol=(0.d0,0.d0)
c ------------------------------------------------------------------
        lam_max=maxval(abs(fif))
        conv_table(0,:,iatom)=lam_max/nproc_t
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          nit=iter_w
          nom=num_omega1_adapt(i_nu)
          nomnu=num_nu2_adapt(i_nu)
          allocate(lambda_dyn(n,n,nom,2,nspin))
          allocate(vrt_dyn_old(n,n,nom,2,nspin))
          allocate(g_om_nu(n,n,2,nom,nspin))
          g_om_nu=(0.d0,0.d0)
c ------ G-interpolation ------------------------------------------
          call g_interp_loc(gexa,gspl_geom,gspl_asy,gasy,g_om_nu,i_nu,
     &                      nom,n)
c --------------------------------------------------------------------
          do s2=1,np
            lambda_dyn=(0.d0,0.d0)
c -------- Get zero order generalized susceptibility -------------------
            call k0_loc(ind_nu,k0_pw,g_om_nu,fif,s2,nom,n,np)
c ----------------------------------------------------------------------
            k_pw=(0.d0,0.d0)
            vrt_dyn_old=(0.d0,0.d0)
            vrt_stat_old=(0.d0,0.d0)
            do it_vrt=1,nit
              k_pw=k0_pw+k_pw
              spin_flips=.false.
              if(it_vrt<=iter_t) spin_flips=.true.
              call vertex_loc_loop(spin_flips,ind_nu,it_vrt,2,s2,
     &                             nom,nomnu,k_pw,fif,vloc(1,1,iatom),
     &                             wloc(1,1,0,iatom),w_tau,lambda_stat,
     &                             lambda_dyn,nspin,n,np,gtau)
c --------------------------------------------------------------------
              conv_1=0.d0
              do ispin=1,nspin
                do i_omega=1,nom
                  do it=1,2
                    do j=1,n
                      do i=1,n
                        cc=lambda_dyn(i,j,i_omega,it,ispin)
     &                    +lambda_stat(i,j,ispin)
     &                    -vrt_dyn_old(i,j,i_omega,it,ispin)
     &                    -vrt_stat_old(i,j,ispin)
                        conv_1=max(conv_1,abs(cc))
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              if(i_nu<=7) then
                conv_table(it_vrt,i_nu,iatom)=
     &              max(conv_1,conv_table(it_vrt,i_nu,iatom))
              endif
              vrt_dyn_old=lambda_dyn
              vrt_stat_old=lambda_stat
              do ispin=1,nspin
                if(it_vrt<nit) then
                  k_pw(:,:,:,:,ispin)=(0.d0,0.d0)
c -------- Get correction to the generalized susceptibility --------
                  call k_loc(ispin,nom,ind_nu,
     &                       lambda_dyn(1,1,1,1,ispin),
     &                       k_pw(1,1,0,1,ispin),
     &                       lambda_stat(1,1,ispin),
     &                       g_om_nu(1,1,1,1,ispin),n)
                else if(it_vrt==nit) then
c -------------- Correction to the Polarizability ---------------------
                  call dp_loc(nom,ind_nu,
     &                        lambda_dyn(1,1,1,1,ispin),
     &                        lambda_stat(1,1,ispin),
     &                        g_om_nu(1,1,1,1,ispin),
     &                        fif(1,1,1,1,ispin),
     &                        pol(1,s2,ind_nu),ispin,
     &                        trn(1,i_nu),1,s2,n,np)
                endif
              enddo  !! over ispin
            enddo  !! over it_vrt
          enddo   !! over s2
          deallocate(lambda_dyn,vrt_dyn_old,g_om_nu)
        enddo   !! over ind_nu
        if(nspin==1.and.irel/=2) pol=2.d0*pol
        ind=iopb_red(iatom)-1
        p_herm=0.d0
        allocate(poli(np,np))
        do ind_nu=1,ndim3_nu
          ph=check_hermit(pol(1,1,ind_nu),np,np)
          p_herm=max(p_herm,ph)
          call pack_hermit(pol(1,1,ind_nu),poli,np,np,np,0.d0,1.d0)
          do ind_q=1,ndim_k_red(me_k+1)
            do j=1,np
              do i=1,np
                p_vertex_nu_q(ind+i,ind+j,ind_nu,ind_q)=poli(i,j)
              enddo
            enddo
          enddo
        enddo
        deallocate(poli)
        if(nproc_t/=1) call dgop(p_herm,1,'max',comm_t)
        if(maswrk) write(iun,*)' VERTEX_P_LOC : IATOM P_HERM = ',iatom,
     &                         p_herm
        deallocate(pol,fif,w_tau,gtau,gexa,gspl_geom,gspl_asy,gasy,
     &             k0_pw,k_pw,lambda_stat,vrt_stat_old)
      enddo  !! over iatom
      deallocate(z_red,trn,g_omega,fifim,g_tau,wloc,vloc)
c ------------------------------------------------------------------
      idum=(iter_w+1)*8*natom
      if(nproc_t/=1) call dgop(conv_table,idum,'  +',comm_t)
      call output_bs_conv_loc(conv_table,iter_w)
      deallocate(conv_table)
      call timel('**** VERTEX_P_LOC finished *********')
      end
