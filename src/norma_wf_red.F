      subroutine norma_wf_red
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: isort,iatom,l,ie,lm,km,i,lget,je,kmj,in,jn,in1,jn1,mj,
     &           li,ib,i_len,ispin,ind_k,k,nbas,j,j0,jj,i0,ii,kc,ka,
     &           iedos,n,isp,npr,ind_b
      real*8 :: zz,deltae,om0,a,emi,ema
      complex*16 :: cc0,cc
      integer, allocatable :: nd(:),n0(:)
      real*8, allocatable :: nrm(:,:,:),gb(:,:),nrl(:),spf(:,:)
      complex*16, allocatable :: s(:,:)
      npr=nproc_t*nproc_b
      allocate(nd(npr))
      allocate(n0(npr))
      allocate(nrm(nbndf,npnt_c,nspin))
      allocate(gb(3,nplwbas))
      nrm=0.d0
      do ispin=1,nspin
        isp=min(ispin,nspin_1)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          kc=k_c_from_a(k)
          if(kc==0) cycle
          nbas=nbask(k)
          do i=1,nbas/nrel
            gb(:,i)=gbs(:,indgb(i,k))
          enddo
          allocate(s(nbas,nbas))
c ------ We form the overlap matrix -----------------------------
          call overlap_lapw(gb,pnt(1,k),nbas,s,indgb(1,k))
          call size_shift_par(n_bnd(k,ispin),npr,nd,n0)
          do ind_b=1,nd(me_t_b+1)
            ib=n0(me_t_b+1)+ind_b
c ------- MT contribution ------------------------------------------
            do iatom=1,natom
              isort=is(iatom)
              do lm=1,nrel*(lmb_red(isort)+1)**2
                if(irel.ne.2) then
                  l=lget(lm)
                  li=l+1
                else if(irel.eq.2) then
                  call getlimj(lm,l,i,mj,li,0)
                endif
                do je=1,ntle(l,isort)
                  if(correlated(je,l,isort)=='N') cycle
                  in1=1
                  if(augm(je,l,isort)=='LAW') in1=2
                  if(augm(je,l,isort)=='AWL') in1=2
                  if(correlated(je,l,isort)=='A') in1=1
                  do jn1=1,in1
                    kmj=io_lem(iatom)-1+indbasa(jn1,je,lm,isort)
                    do ie=1,ntle(l,isort)
                      if(correlated(ie,l,isort)=='N') cycle
                      in=1
                      if(augm(ie,l,isort)=='LAW') in=2
                      if(augm(ie,l,isort)=='AWL') in=2
                      if(correlated(ie,l,isort)=='A') in=1
                      do jn=1,in
                        km=io_lem(iatom)-1+indbasa(jn,ie,lm,isort)
                        zz=conjg(z_bnd(km,ib,ind_k,ispin))
     &		  	              *z_bnd(kmj,ib,ind_k,ispin)
                        nrm(ib,kc,ispin)=nrm(ib,kc,ispin)
     &	                    +zz*ffsmt(jn,jn1,ie,je,li,isort,isp)
                      enddo
                    enddo
                  enddo
                enddo
              enddo    !! over lm
            enddo  !! over iatom
c ------- Interstitial contribution ------------------------------
            do j=1,nbask_red(k)/nrel  !! in reduced basis for k
              j0=indgb_red(j,k)  !! in GBS-list
              jj=iplf_bk(j0,k)    !! in full basis for k
              do i=1,nbask_red(k)/nrel  !! in reduced basis for k
                i0=indgb_red(i,k)  !! in GBS-list
                ii=iplf_bk(i0,k)    !! in full basis for k
                zz=conjg(ev_bnd(ii,ib,ind_k,ispin))
     &		  	      *ev_bnd(jj,ib,ind_k,ispin)
                nrm(ib,kc,ispin)=nrm(ib,kc,ispin)
     &	                      +zz*s(ii,jj)
              enddo
            enddo
          enddo   !! over ind_b
          deallocate(s)
        enddo   !! over ind_k
      enddo   !! over ispin
      if(nproc_k/=1) call dgop(nrm,nbndf*npnt_c*nspin,'  +',comm_k)
      if(npr/=1) call dgop(nrm,nbndf*npnt_c*nspin,'  +',comm_t_b)
      allocate(nrl(nc_line))
      i_len=len_trim(allfile)
      if(maswrk) then
        do ispin=1,nspin
          if(ispin==1) open(3,file=allfile(1:i_len)//'_NORMA_WF_RED_DN')
          if(ispin==2) open(3,file=allfile(1:i_len)//'_NORMA_WF_RED_UP')
          do ib=1,nbndf
            do kc=1,nc_line
              k=kline_in_npnt_c(kc)
              nrl(kc)=nrm(ib,k,ispin)
            enddo
            write(3,'(i4,5(1x,e12.5))')ib,(nrl(k),k=1,min(nc_line,5))
          enddo
          close(3)
        enddo
      endif
      deallocate(nrl,gb)
! --------- X-Spectral function in RED basis and NPNT_C mesh -----------
      emi=1000.d0
      ema=-1000.d0
      do ispin=1,nspin
        do k=1,npnt_c
          ka=k_a_from_c(k)
          do ib=1,n_low_bnd(ka,nspin)
            i=ind_bands_bnd(ib,ka,ispin)
            emi=min(emi,e_bnd(i,ka,ispin)-chem_pot)
            ema=max(ema,e_bnd(i,ka,ispin)-chem_pot)
          enddo
        enddo
      enddo
      emi=1.3d0*emi
      ema=1.3d0*ema
      n=(ema-emi)/(emaxdos-emindos)*ndos
      deltae=(ema-emi)/n
      allocate(spf(0:n,2))
      spf=0.d0
      do ispin=1,nspin
        do ind_k=1,ndim_k_red(me_k+1)
          k=n_mpi_k_red(me_k+1)+ind_k
          ka=k_a_from_c(k)
          call size_shift_par(n_bnd(ka,ispin),npr,nd,n0)
          do ind_b=1,nd(me_t_b+1)
            ib=n0(me_t_b+1)+ind_b
            cc0=dcmplx(chem_pot-e_bnd(ib,ka,ispin),e_small*10.d0)
            do iedos=0,n
              om0=emi+iedos*deltae
              cc=cc0+om0
              a=imag((1.d0,0.d0)/cc)*wgt_c(k)/pi
              spf(iedos,2)=spf(iedos,2)-a*nrm(ib,k,ispin)
              spf(iedos,1)=spf(iedos,1)-a
            enddo
          enddo
        enddo
      enddo
      if(nproc_k/=1) call dgop(spf,2*(n+1),'  +',comm_k)
      if(npr/=1) call dgop(spf,2*(n+1),'  +',comm_t_b)
      spf=2.d0*spf/nspin/nrel
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'X_FULL_RED.dos')
        do iedos=0,n
          om0=emi+iedos*deltae
          write(3,'(3(f11.5,2x))')om0*evolt/2,spf(iedos,1),spf(iedos,2)
        enddo
        close(3)
      endif
      deallocate(nrm,spf,n0,nd)
      call timel('**** NORMA_WF_RED finished *********')
      end
