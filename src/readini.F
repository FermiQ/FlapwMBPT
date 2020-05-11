	  subroutine readini
	  use atom_mod
	  use etot_mod
	  use heg_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
	  implicit none
#ifdef MPI
	  include 'mpif.h'
#endif
	  character*1 :: chh
	  character*20 :: txt
	  character*80 :: cv
	  integer :: iatom,i,it,isort,lb,lpb,l,ll,n,ntl(0:20),i1,ii
	  real*8 :: aa,bb,cc,alf,bet,gam,sgam,calf,cbet,cgam,salf
	  if(maswrk) then
	    open(1,file='ini',status='old')
	    read(1,'(a1)')chh
	    read(1,*)txt,txt,iter_dft,txt,iter_hf,txt,iter_gw,txt,iter_qp
	    read(1,*)txt,iter_psi,txt,iter_bsp
	    read(1,*)txt,restart_begin,txt,restart_end
	    read(1,*)txt,admix,txt,adspin,txt,adm_gw,txt,acc_it_gw
	    read(1,*)txt,iexch,txt,scal_spin,txt,psi_fncl_use
	    read(1,*)txt,irel,txt,clight,txt,rel_interst
	    nrel=1
	    if(irel.eq.2) nrel=2
	    read(1,*)txt,temperature
	    read(1,*)txt,txt,mpi_pref
        read(1,*)txt,nproc_t,txt,nproc_k,txt,nproc_b
	    read(1,'(a1)')chh
	    read(1,*)txt,allfile
        read(1,*)txt,tmpfile
	    read(1,*)txt,txt,symgen
	    read(1,*)txt,txt,par,txt,nsort,txt,istruc
	  endif
	  if(goparr) then
#ifdef MPI
	    call brdcst(msgint,iter_dft,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iter_psi,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iter_bsp,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iter_hf,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,iter_gw,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iter_qp,4,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,admix,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,adspin,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,clight,8,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iexch,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,irel,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nrel,4,master,MPI_COMM_WORLD)
	    call brdcstc(allfile,72,master,MPI_COMM_WORLD)
        call brdcstc(tmpfile,72,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nproc_t,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nproc_k,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,rel_interst,4,master,MPI_COMM_WORLD)
	    call brdcstc(symgen,80,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,par,8,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nsort,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,istruc,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nproc_b,4,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,scal_spin,8,master,MPI_COMM_WORLD)
	    call brdcst(msgint,restart_begin,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,restart_end,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,psi_fncl_use,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,adm_gw,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,acc_it_gw,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,temperature,8,master,MPI_COMM_WORLD)
        call brdcstc(mpi_pref,3,master,MPI_COMM_WORLD)
#endif	  
	  endif
	  mbpt_gw=0
	  if(iter_gw<0) then
        mbpt_gw=iabs(iter_gw)
	    iter_gw=1
	  endif
	  allocate(tau(3,nsort))
	  allocate(txtel(nsort))
	  allocate(z(nsort))
	  allocate(z_dop(nsort))
	  allocate(zcor(nsort))
	  allocate(magn_shift(nsort))
	  allocate(smt(nsort))
	  allocate(h(nsort))
	  allocate(nrad(nsort))
	  allocate(lmb(nsort))
	  allocate(lmb_red(nsort))
	  allocate(lmpb(nsort))
      allocate(lmbs(nsort))
	  if(maswrk) then
	    read(1,*)txt,b_a,txt,c_a
	    read(1,*)txt,cv
	    i1=index(cv,' ')-1
        call decode_vector(cv,1,i1,rbas(1,1))
        read(1,*)txt,cv
        i1=index(cv,' ')-1
        call decode_vector(cv,1,i1,rbas(1,2))
        read(1,*)txt,cv
        i1=index(cv,' ')-1
        call decode_vector(cv,1,i1,rbas(1,3))
        aa=dot_product(rbas(:,3),rbas(:,3))
        if(aa<1.d-12) then
          if(iabs(istruc)==12) then
c --------- Use a,b,c and 3 angles for simple monoclinic -----------
            aa=rbas(1,1)
            bb=rbas(2,1)
            cc=rbas(3,1)
            alf=rbas(1,2)*pi/180.d0
            salf=sin(alf)
            calf=cos(alf)
            rbas=0.d0
            rbas(1,1)=aa
            rbas(2,2)=bb
            rbas(2,3)=cc*calf
            rbas(3,3)=cc*salf
          else if(iabs(istruc)==15) then
c --------- Use a,b,c and 3 angles for triclinic solid -------------
            aa=rbas(1,1)
            bb=rbas(2,1)
            cc=rbas(3,1)
            alf=rbas(1,2)*pi/180.d0
            bet=rbas(2,2)*pi/180.d0
            gam=rbas(3,2)*pi/180.d0
            sgam=sin(gam)
            calf=cos(alf)
            cbet=cos(bet)
            cgam=cos(gam)
            rbas=0.d0
            rbas(1,1)=aa
            rbas(1,2)=bb*cgam
            rbas(2,2)=bb*sgam
            rbas(1,3)=cc*cbet
            rbas(2,3)=cc*(calf-cbet*cgam)/sgam
            rbas(3,3)=cc*sqrt(sgam**2-calf**2-cbet**2
     &                     +2.d0*calf*cbet*cgam)/sgam
          endif
        endif
	    do iatom=1,nsort
          if(iatom==1) read(1,*)txt,cv
          if(iatom/=1) read(1,*)cv
          i1=index(cv,' ')-1
          call decode_vector(cv,1,i1,tau(1,iatom))
	    enddo
	    read(1,*)txt,txt,txt,txt,mdiv
	    read(1,*)txt,nrdiv
	    read(1,*)txt,nrdiv_red
	    read(1,*)txt,txt,cut_gw_ratio,txt,cut_vrt_ratio
	    read(1,*)txt,eps_pb,txt,eps_pb_vrt
	    read(1,*)txt,txt,nbndf,txt,nbndf_bnd_0
	    read(1,*)txt,txt,n_k_div
	    read(1,*)txt,txt,emindos,txt,emaxdos,txt,ndos
	    read(1,*)txt,e_small,txt,e_small_bos
	    read(1,*)txt,txt,(ndiv(i),i=1,3),txt,k_line
	    read(1,*)txt,(ndiv_c(i),i=1,3)
        read(1,*)txt,k_integral,txt,(n_k_int(i),i=1,3)
	    read(1,*)txt,txt,v_v0
	    read(1,*)txt,txt,b_extval,txt,iter_h_ext
c -------- Transform to internal magnetic field --------------
	    read(1,*)txt,(b_ext(i),i=1,3)
	    b_ext=b_ext/sqrt(dot_product(b_ext,b_ext))
	    read(1,*)txt,txt,eps_coul
	    read(1,*)txt,txt,w_sc_gw,txt,w_sc_qp,txt,w_sc_psi,txt,w_sc_bsp
        read(1,*)txt,txt,alpha_hbr,txt,omega_hse,txt,hybrid_type
	    read(1,*)txt,txt,psi_p,txt,bse_kernel_p
	    read(1,*)txt,iter_ladder_p,txt,nu_w_stat
	    read(1,*)txt,txt,psi_sig,txt,psi2_sig
	    read(1,*)txt,iter_sigma_gwg
	    read(1,*)txt,txt,iter_ladder_xi,txt,q_suscept
	    read(1,*)txt,chi_cmp,txt,vrt_x_appr
	    read(1,*)txt,nrax_chi,txt,freq_chi
	    read(1,*)txt,theta_bse
        read(1,*)txt,txt,opt_mode,txt,opt_loc_fld
	    read(1,*)txt,txt,txt,n_tau,txt,exp_tau_gw
	    read(1,*)txt,txt,txt,n_omega_exa,txt,n_omega_geom,txt,
     &           n_omega_asy
	    read(1,*)txt,omega_geom,txt,omega_max
	    omega_last_step=1
	    read(1,*)txt,txt,txt,n_nu_exa,txt,n_nu_geom,txt,n_nu_asy
	    read(1,*)txt,nu_geom,txt,nu_max
	    nu_last_step=1
	    read(1,*)txt,txt,rs_el_gas,txt,k_inf,txt,n_k_heg
        n_bnd_out=1
        inf_bnd_out(:,1)=(/1,1,1/)
	    read(1,'(a1)')chh
	    maxb=0
	    maxb_red=0
	    maxpb=0
	    maxntle=0
	    maxs=0
	    it=0
	    do isort=1,nsort
		  read(1,'(a1)')chh
		  read(1,'(a1)')chh
		  read(1,*)txt,lb,txt,lpb
		  read(1,'(a1)')chh
		  read(1,'(a1)')chh
		  maxb=max(maxb,lb)
		  maxpb=max(maxpb,lpb)
		  maxs=max(maxs,lpb+8)
		  read(1,*)txt,(ntl(l),l=0,lb)
		  read(1,'(a1)')chh
		  it=it+7
		  do l=0,lb
		    do n=1,ntl(l)
			  it=it+1
			  read(1,'(a1)')chh
		    enddo
		    maxntle=max(maxntle,abs(ntl(l)))
		  enddo
	    enddo
	  endif
	  if(goparr) then
#ifdef MPI
	    call brdcst(msgint,n_omega_exa,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_omega_geom,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_omega_asy,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,omega_last_step,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_nu_exa,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_nu_geom,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nu_last_step,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_nu_asy,4,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,omega_geom,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,omega_max,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,exp_tau_gw,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,nu_geom,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,nu_max,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,cut_gw_ratio,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,cut_vrt_ratio,8,master,MPI_COMM_WORLD)
	    call brdcst(msgint,maxb,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,maxpb,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,maxs,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,maxntle,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_k_div,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iter_ladder_xi,4,master,MPI_COMM_WORLD)
	    call brdcstc(k_line,3,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_gamma,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iter_h_ext,4,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,eps_pb_vrt,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,e_small,8,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,e_small_bos,8,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_tau,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iter_ladder_p,4,master,MPI_COMM_WORLD)
	    call brdcstc(vrt_x_appr,3,master,MPI_COMM_WORLD)
	    call brdcstc(w_sc_gw,3,master,MPI_COMM_WORLD)
	    call brdcstc(w_sc_qp,3,master,MPI_COMM_WORLD)
	    call brdcstc(w_sc_psi,3,master,MPI_COMM_WORLD)
	    call brdcstc(w_sc_bsp,3,master,MPI_COMM_WORLD)
	    call brdcstc(psi_sig,5,master,MPI_COMM_WORLD)
	    call brdcstc(psi2_sig,5,master,MPI_COMM_WORLD)
	    call brdcstc(psi_p,5,master,MPI_COMM_WORLD)
	    call brdcstc(bse_kernel_p,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,iter_sigma_gwg,4,master,MPI_COMM_WORLD)
	    call brdcstc(chi_cmp,3,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nrax_chi,4,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,freq_chi,8,master,MPI_COMM_WORLD)
	    call brdcst(msgint,ndiv,4*3,master,MPI_COMM_WORLD)
	    call brdcst(msgint,mdiv,4*3,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nrdiv,4*3,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nrdiv_red,4*3,master,MPI_COMM_WORLD)
	    call brdcst(msgint,ndiv_c,4*3,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,k_inf,8,master,MPI_COMM_WORLD)
	    call brdcst(msgint,n_k_heg,4,master,MPI_COMM_WORLD)
	    call brdcstc(theta_bse,3,master,MPI_COMM_WORLD)
	    call brdcst(msgint,q_suscept,4,master,MPI_COMM_WORLD)
	    call brdcst(msgint,nu_w_stat,4,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,alpha_hbr,24,master,MPI_COMM_WORLD)
	    call brdcstc(opt_mode,6,master,MPI_COMM_WORLD)
	    call brdcst(msgint,opt_loc_fld,4,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,eps_coul,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_bnd_out,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,inf_bnd_out,72,master,MPI_COMM_WORLD)
        call brdcstc(k_integral,2,master,MPI_COMM_WORLD)
        call brdcst(msgint,n_k_int,12,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,rs_el_gas,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,b_a,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,c_a,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,rbas,8*9,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,tau,8*3*nsort,master,MPI_COMM_WORLD)
        call brdcst(msgint,nbndf,4,master,MPI_COMM_WORLD)
        call brdcst(msgint,nbndf_bnd_0,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,ndos,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,emindos,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,emaxdos,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,v_v0,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,b_extval,8,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,b_ext,8*3,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,eps_pb,8,master,MPI_COMM_WORLD)
        call brdcst(msgint,hybrid_type,4,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,omega_hse,8,master,MPI_COMM_WORLD)
#endif	  
	  endif
c ----- TAU in Cartesian coordinates if necessary ------------------
      if(istruc.lt.0) then
        istruc=-istruc
        do iatom=1,nsort
          aa=tau(1,iatom)
          bb=tau(2,iatom)
          cc=tau(3,iatom)
          do i=1,3
            tau(i,iatom)=aa*rbas(i,1)+bb*rbas(i,2)+cc*rbas(i,3)
          enddo
        enddo
      endif
c ------------------------------------------------------------------
	  allocate(konfig(0:maxb,nsort))
	  allocate(atoc(maxntle,0:maxb,nsort))
	  allocate(idmd(maxntle,0:maxb,nsort))
	  allocate(ntle(0:maxb,nsort))
	  allocate(augm(maxntle,0:maxb,nsort))
	  allocate(correlated(maxntle,0:maxb,nsort))
	  allocate(ptnl(maxntle,0:maxb,nsort,2))
	  allocate(lim_pb_mt(0:maxpb,nsort))
	  allocate(lim_pb_mt_red(0:maxpb,nsort))
	  h=0.d0
	  nrad=0
	  z=0.d0
	  z_dop=0.d0
	  atoc=0.d0
	  if(maswrk) then
	    do i=1,it
		  backspace 1
	    enddo
	    lmb_red=-1
	    do isort=1,nsort
          read(1,*)txt,txtel(isort),txt,z(isort),txt,magn_shift(isort),
     &             txt,z_dop(isort)
		  read(1,*)txt,smt(isort),txt,h(isort),txt,nrad(isort)
		  read(1,*)txt,lmb(isort),txt,lmpb(isort)
		  lmbs(isort)=lmpb(isort)+8
		  read(1,*)txt,(lim_pb_mt(i,isort),i=0,lmpb(isort))
		  read(1,*)txt,(lim_pb_mt_red(i,isort),i=0,lmpb(isort))
		  read(1,*)txt,(ntle(l,isort),l=0,lmb(isort))
		  read(1,'(a1)')chh
c ------- Reading the type of the functions for L <= L_b ------
		  do l=0,lmb(isort)
		    konfig(l,isort)=10
		    do n=1,ntle(l,isort)
			  read(1,*)ll,augm(n,l,isort),atoc(n,l,isort),
     &                 ptnl(n,l,isort,1),correlated(n,l,isort),
     &                 idmd(n,l,isort)
	   	      ptnl(n,l,isort,2)=ptnl(n,l,isort,1)
			  konfig(l,isort)=min(konfig(l,isort),
     &                            int(ptnl(n,l,isort,1)))
		      if(correlated(n,l,isort)/='N')
     &                        lmb_red(isort)=max(l,lmb_red(isort))
		    enddo
		  enddo
	    enddo
	    close(1)
	    maxb_red=maxval(lmb_red)
	  endif
	  if(goparr) then
#ifdef MPI
	    call brdcstc(txtel,4*nsort,master,MPI_COMM_WORLD)
	    call dgop(z,nsort,'  +',MPI_COMM_WORLD)
	    call dgop(z_dop,nsort,'  +',MPI_COMM_WORLD)
	    call brdcst(msgdbl,magn_shift,8*nsort,master,MPI_COMM_WORLD)
	    call brdcst(msgdbl,smt,8*nsort,master,MPI_COMM_WORLD)
	    call dgop(h,nsort,'  +',MPI_COMM_WORLD)
	    call igop(nrad,nsort,'  +',MPI_COMM_WORLD)
	    call dgop(atoc,nsort*(maxb+1)*maxntle,'  +',MPI_COMM_WORLD)
	    call brdcstc(correlated,nsort*(maxb+1)*maxntle,master,
     &               MPI_COMM_WORLD)
	    call brdcst(msgint,idmd,4*nsort*(maxb+1)*maxntle,master,
     &              MPI_COMM_WORLD)
	    call brdcst(msgdbl,ptnl,16*nsort*(maxb+1)*maxntle,master,
     &              MPI_COMM_WORLD)
	    call brdcstc(augm,3*nsort*(maxb+1)*maxntle,master,
     &               MPI_COMM_WORLD)
	    call brdcst(msgint,lmb,4*nsort,master,MPI_COMM_WORLD)
        call brdcst(msgint,lmbs,4*nsort,master,MPI_COMM_WORLD)
	    call brdcst(msgint,lmb_red,4*nsort,master,MPI_COMM_WORLD)
	    call brdcst(msgint,lmpb,4*nsort,master,MPI_COMM_WORLD)
	    call brdcst(msgint,lim_pb_mt,4*nsort*(maxpb+1),master,
     &              MPI_COMM_WORLD)
	    call brdcst(msgint,lim_pb_mt_red,4*nsort*(maxpb+1),master,
     &              MPI_COMM_WORLD)
	    call brdcst(msgint,konfig,4*nsort*(maxb+1),master,
     &              MPI_COMM_WORLD)
        call brdcst(msgint,ntle,4*nsort*(maxb+1),master,MPI_COMM_WORLD)
	    call brdcst(msgint,maxb_red,4,master,MPI_COMM_WORLD)
#endif	  
	  endif
	  z=z+z_dop
c ------------------------------------------------------------------
      if(maswrk) then
        i1=len_trim(allfile)
        do ii=i1+1,72
          dosfile(ii:ii)=' '
          outfile(ii:ii)=' '
          rstfile(ii:ii)=' '
        enddo
        dosfile(1:i1+4)=allfile(1:i1)//'.dos'
        outfile(1:i1+4)=allfile(1:i1)//'.out'
        rstfile(1:i1+4)=allfile(1:i1)//'.rst'
      endif
c     Output file
      iun=2
      if(maswrk) then
        open(iun,file=outfile)
        call goto_end(iun)
      endif
	  end
