      subroutine afteread
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      use units_mod
      implicit none
      integer :: iatom,i,isort,iz,l,nl,nlm,n,ind,li,ie,ii,jj,ii2,ndim,j,
     &           lm,m,ndimb,iis
      real*8 :: am(103),aa,ro_el_0
      DATA Am/  1.0079D0,  4.0026D0,  6.941D0,  9.01218D0,  10.81D0,
     &        12.011D0,   14.0067D0,  15.999D0, 18.998403D0,   20.17D0,
     &        22.98977D0, 24.305D0,  26.98154D0, 28.085D0,   30.97376D0,
     &        32.064D0,  35.453D0,  39.948D0, 39.098D0,    40.08D0,
     &        44.9559D0,  47.90D0,   50.9415D0, 51.996D0,    54.938D0,
     &        55.84D0,   58.9332D0,   58.70D0,  63.546D0,    65.38D0,
     &        69.72D0,   72.59D0,   74.9216D0, 78.96D0,    79.904D0,
     &        83.80D0,   85.4678D0,   87.62D0,  88.9059D0,    91.22D0,
     &        92.9064D0,   95.94D0,   98.91D0, 101.07D0,   102.9055D0,
     &       106.40D0,  107.868D0,  112.41D0, 114.82D0,   118.69D0,
     &       121.75D0,  127.60D0,  126.9045D0, 131.30D0,   132.9054D0,
     &       137.33D0,  138.9055D0,  140.12D0, 140.9077D0,   144.24D0,
     &       145.D0,    150.35D0,  151.96D0, 157.25D0,   158.9254D0,
     &       162.50D0,  164.9304D0,  167.26D0, 168.9342D0,   173.04D0,
     &       174.967D0,  178.49D0,  180.9479D0, 183.85D0,   186.207D0,
     &       190.20D0,  192.22D0,  195.09D0, 196.9665D0,   200.59D0,
     &       204.37D0,  207.19D0,  208.9804D0, 209.D0,     210.D0,
     &       222.D0,    223.D0,    226.0254D0,   227.0278D0, 232.0381D0,
     &       231.0359D0, 238.029D0,  237.0482D0, 244.D0,     243.D0,
     &       247.D0,    247.D0,    251.D0,   252.D0,     257.D0,
     &       258.D0,    259.D0,    260.D0/
c     *******************************************************
      allocate(lm_general(nrel*maxb+1,2*maxb+nrel))
      lm_general=0
      lm=0
      do li=1,nrel*maxb+1
        if(irel.ne.2) ndim=li+li-1
        if(irel.eq.2) ndim=((li+1)/2)*2
        do m=1,ndim
          lm=lm+1
          lm_general(li,m)=lm
        enddo
      enddo
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      molar_mass=0.d0
      do isort=1,nsort
         iz=idnint(z(isort))
         if(iz>0) then
           amass(isort)=am(iz)*911.43957d0
           molar_mass=molar_mass+am(iz)*nhsort(isort)
         else
           amass(isort)=0.d0
         endif
      enddo
c ------- Core information --------------------------------------------
      ncormax=0
      ncormmax=0
      ncor_tot=0
      do isort=1,nsort
        zcor(isort)=0.d0
        ncor(isort)=0
        ncorm(isort)=0
        lc_max(isort)=-1
        do l=0,lmb(isort)
          nl=1
          if(irel.eq.2.and.l.ne.0) nl=2
          nlm=nrel*(l+l+1)
          do n=l+1,konfig(l,isort)-1
            zcor(isort)=zcor(isort)+2.d0*(l+l+1.d0)
            ncor(isort)=ncor(isort)+nl
            ncorm(isort)=ncorm(isort)+nlm
            lc_max(isort)=max(lc_max(isort),l)
          enddo
        enddo
        ncormax=max(ncormax,ncor(isort))
        ncormmax=max(ncormmax,ncorm(isort))
        ncor_tot=ncor_tot+ncorm(isort)*nhsort(isort)
      enddo
      maxl_cor=maxval(lc_max)
      ind=1
      do iatom=1,natom
        isort=is(iatom)
        io_cor(iatom)=ind
        ind=ind+ncorm(isort)
      enddo
c
c                 NELEC - The number of valence electrons
c
      nelec=0.d0
      do isort=1,nsort
        aa=0.d0
        do li=1,nrel*lmb(isort)+1
          if(irel.le.1) then
            l=li-1
          else
            l=li/2
          endif
          do ie=1,ntle(l,isort)
            if(l==0) aa=aa+atoc(ie,l,isort)
            if(l>0) aa=aa+atoc(ie,l,isort)/nrel
          enddo
        enddo
        nelec=nelec+(aa+z_dop(isort))*nhsort(isort)
      enddo
c ------- Valence information --------------------------------------------
      maxel=0
      maxlfun=0
      nfun=0
      nfun0=0
      nlb=0
      nlb2=0
      nlbs=0
      nlb20=0
      do isort=1,nsort
        lfun(isort)=0
        lfunm(isort)=0
        do li=1,nrel*lmb(isort)+1
          if(irel.ne.2) then
            l=li-1
            nlm=l+l+1
          else if(irel.eq.2) then
            l=li/2
            nlm=((li+1)/2)*2
          endif
          do ie=1,ntle(l,isort)
            lfun(isort)=lfun(isort)+1
            lfunm(isort)=lfunm(isort)+nlm
            if(augm(ie,l,isort)=='LOC'.or.augm(ie,l,isort)=='APW') then
              indfun0(1,ie,li,isort)=lfun(isort)
            else if(augm(ie,l,isort)=='LAW'.or.
     &              augm(ie,l,isort)=='AWL') then
              indfun0(1,ie,li,isort)=lfun(isort)
              lfun(isort)=lfun(isort)+1
              lfunm(isort)=lfunm(isort)+nlm
              indfun0(2,ie,li,isort)=lfun(isort)
            endif
          enddo
        enddo
        maxel=max(maxel,lfunm(isort))
        maxlfun=max(maxlfun,lfun(isort))
        nfun=nfun+lfunm(isort)*nhsort(isort)
        nfun0=nfun0+lfunm(isort)
        nlb=nlb+(lmb(isort)+1)**2*nhsort(isort)
        nlb2=nlb2+(lmpb(isort)+1)**2*nhsort(isort)
        nlbs=nlbs+(lmbs(isort)+1)**2*nhsort(isort)
        nlb20=nlb20+(lmpb(isort)+1)**2
      enddo
      allocate(iat_tel(nfun))
      allocate(iso_tel(nfun0))
c ---------------------------------------------------------
      ii=1
      jj=1
      ii2=1
      iis=1
      do iatom=1,natom
        isort=is(iatom)
        iob(iatom)=ii
        iob2(iatom)=ii2
        iobs(iatom)=iis
        io_lem(iatom)=jj
        ii=ii+(lmb(isort)+1)**2
        ii2=ii2+(lmpb(isort)+1)**2
        iis=iis+(lmbs(isort)+1)**2
        iat_tel(jj:jj+lfunm(isort)-1)=iatom
        jj=jj+lfunm(isort)
      enddo
c ----------------------------------------------------
      jj=1
      ii2=1
      do isort=1,nsort
        iob20(isort)=ii2
        io_lem0(isort)=jj
        ii2=ii2+(lmpb(isort)+1)**2
        iso_tel(jj:jj+lfunm(isort)-1)=isort
        jj=jj+lfunm(isort)
      enddo
c ----------------------------------------------------
      nlcloc=0
      max_bs=0
      do iatom=1,natom
        isort=is(iatom)
        nlctot(isort)=0
        ndim=nrel*lmb(isort)+1
        do li=1,ndim
          if(irel.ne.2) l=li-1
          if(irel.eq.2) l=li/2
          do ie=1,ntle(l,isort)
            if(irel.ne.2) nl=l+l+1
            if(irel.eq.2) nl=((li+1)/2)*2
            if(augm(ie,l,isort)=='LOC'.or.augm(ie,l,isort)=='AWL') then
              nlcloc=nlcloc+nl
              nlctot(isort)=nlctot(isort)+nl
            endif
          enddo
        enddo
        max_bs=max_bs+2*nrel*(lmb(isort)+1)**2
      enddo
      i=0
      do iatom=1,natom
        isort=is(iatom)
        lm=0
        do li=1,nrel*lmb(isort)+1
          if(irel.ne.2) l=li-1
          if(irel.eq.2) l=li/2
          if(irel.ne.2) ndim=li+li-1
          if(irel.eq.2) ndim=((li+1)/2)*2
          do m=1,ndim
            lm=lm+1
            do ie=1,ntle(l,isort)
              if(augm(ie,l,isort)=='LOC'.or.
     &           augm(ie,l,isort)=='AWL') then
                i=i+1
                indbas(ie,lm,iatom)=i
              endif
            enddo
          enddo
        enddo
      enddo
c -------  Core information ----------------------------------
      maxnrad=0
      maxmtcor=0
      do isort=1,nsort
        maxmtcor=maxmtcor+(nrad(isort)+1)*ncor(isort)*nspin_0
        maxnrad=max(maxnrad,nrad(isort))
      enddo
c ------ Here are things for the new augmentation procedure ----------
      allocate(n_e_apw(0:maxb,nsort))
      n_e_apw=0
      do isort=1,nsort
        do l=0,lmb(isort)
          do ie=1,ntle(l,isort)
            if(augm(ie,l,isort)=='APW') n_e_apw(l,isort)=ie
            if(augm(ie,l,isort)=='LAW') n_e_apw(l,isort)=ie
            if(augm(ie,l,isort)=='AWL') n_e_apw(l,isort)=ie
          enddo   !! over ie
        enddo  !! over l
      enddo   !! over isort
c -----------------------------------------------------------     
      allocate(ind_y_apw(2,limlb,natom))
      ind=0
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        do ie=1,2
          do lm=1,ndimb
            ind=ind+1
            ind_y_apw(ie,lm,iatom)=ind
          enddo
        enddo
      enddo
c -------- Number of semicore bands --------------------------------
      n_semic=0
      do isort=1,nsort
        do l=0,lmb(isort)
          ii=n_e_apw(l,isort)  !! augmentation energy
          jj=int(ptnl(ii,l,isort,1))  !! corresp. princ.quant.number
          do ie=1,ntle(l,isort)
            if(augm(ie,l,isort)=='LOC') then
              j=int(ptnl(ie,l,isort,1))
              if(j<jj) n_semic=n_semic+nrel*(l+l+1)*nhsort(isort)
            endif
          enddo
       enddo
      enddo
      iheg=.false.
      if(nelec<1.d-10) then
!     &&&&& Here for electron gas only &&&&&&&&&&&&&&&&&&&&&&&&&&
        ro_el_0=3.d0/4.d0/pi/rs_el_gas**3
        iheg=.true.
        nelec=1.d0
        par=par*(nelec/ro_el_0/amega)**(1.d0/3.d0)
        amega=nelec/ro_el_0
      endif
      if(irel/=2) ro_el_gas=(nelec-n_semic)/amega
      if(irel==2) ro_el_gas=(nelec-n_semic)/amega
      end
