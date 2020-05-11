      subroutine transport_boltz
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: k,i,i_len,ispin,j,nt,ib,it,n,ig,nd(3),ndd,np
      real*8 :: a(3,3),e2t,t,dt,tmi,tma,nee,beta,dum,nel_b,de2,kel,zt1,
     &          fermi_dirac_f,de,f,df,p1(3,3),p2(3,3),per(3,3),zt,s2,s,
     &          e_don,s_don
      real*8, allocatable :: eb(:,:,:),deb(:,:,:,:),sig0(:,:,:,:,:),
     &                       sig(:,:,:),nu(:,:,:),seeb(:,:,:),kp(:,:),
     &                       wg(:),ebnd(:,:,:),chemp(:),res(:,:,:),
     &                       kap0(:,:,:),kap_el(:,:,:)
      complex*16, allocatable :: aa(:,:)
      n_b_t=14
      t_b_min=300.d0
      t_b_max=1000.d0
      tau_relax=3.6d-15
      e2t=tau_relax/2.418884d-17
      tmi=t_b_min
      tma=t_b_max
      nt=n_b_t
      dt=(tma-tmi)/nt
      call read_boltzman(0,n,nel_b,dum,e_don,s_don)
      allocate(ebnd(n,npnt,nspin))
      call read_boltzman(1,n,nel_b,ebnd,e_don,s_don)
      do ispin=1,nspin
        do k=1,npnt
          do i=4,n
            ebnd(i,k,ispin)=ebnd(i,k,ispin)-0.3243d0/evolt*2.d0
          enddo
        enddo
      enddo
c -------------- K - mesh generation -------------------------------
      nd=ndiv*4
      ndd=nd(1)*nd(2)*nd(3)
      allocate(kp(3,ndd))
      allocate(wg(ndd))
      call k_points_gener(nd,ndd,np,kp,wg)
c ------------------------------------------------------------------
      allocate(aa(nstar_r,n*nspin))
      call bnd_interp_prep(n,ebnd,aa)
      allocate(eb(n,np,nspin))
      allocate(deb(n,np,nspin,3))
      call bnd_interp(aa,np,kp,n,eb,deb)
      deallocate(aa)
c ------------------------------------------------------------------
      allocate(sig0(3,3,n,np,nspin))
      do ispin=1,nspin
        do k=1,np
          do ib=1,n
            do j=1,3
              do i=1,3
                sig0(i,j,ib,k,ispin)=e2t*deb(ib,k,ispin,i)
     &                                  *deb(ib,k,ispin,j)
              enddo
            enddo
          enddo
        enddo
      enddo
      nee=nel_b+s_don
      allocate(sig(3,3,0:nt))
      allocate(nu(3,3,0:nt))
      allocate(kap0(3,3,0:nt))
      allocate(chemp(0:nt))
      sig=0.d0
      nu=0.d0
      kap0=0.d0
      do it=0,nt
        t=tmi+dt*it
        beta=1.d0/t/boltzman
        call search_mu_trnsp(chem_pot,nee,chemp(it),eb,n,beta,np,wg,
     &                       e_don,s_don)
c        chemp(it)=chemp(0)
        do ispin=1,nspin
          do k=1,np
            do ib=1,n
              de=eb(ib,k,ispin)-chemp(it)
              f=fermi_dirac_f(de*beta)
              df=f*(1.d0-f)*wg(k)
              de2=de**2*df
              de=de*df
              do j=1,3
                do i=1,3
                  sig(i,j,it)=sig(i,j,it)+sig0(i,j,ib,k,ispin)*df
                  nu(i,j,it)=nu(i,j,it)+sig0(i,j,ib,k,ispin)*de
                  kap0(i,j,it)=kap0(i,j,it)+sig0(i,j,ib,k,ispin)*de2
                enddo
              enddo
            enddo
          enddo  !!! over k
        enddo  !!! over ispin
        sig(:,:,it)=sig(:,:,it)/t
        nu(:,:,it)=-nu(:,:,it)/t/t
        kap0(:,:,it)=kap0(:,:,it)/t/t
      enddo
      sig=sig*2.d0/nspin/nrel/amega !! a.u. Rydberg
      nu=nu*2.d0/nspin/nrel/amega  !! a.u. Rydberg
      kap0=kap0*2.d0/nspin/nrel/amega  !! a.u. Rydberg
c -------- Symmetrization ------------------------------------------
      do i=0,nt
        per=0.d0
        do ig=1,ngroup
          do j=1,3
            call rotate_r(sig(j,:,i),p1(j,:),u(2,ig),1)
          enddo
          do j=1,3
            call rotate_r(p1(1,j),p2(1,j),u(2,ig),1)
          enddo
          per=per+p2/ngroup
        enddo
        sig(:,:,i)=per
        per=0.d0
        do ig=1,ngroup
          do j=1,3
            call rotate_r(nu(j,:,i),p1(j,:),u(2,ig),1)
          enddo
          do j=1,3
            call rotate_r(p1(1,j),p2(1,j),u(2,ig),1)
          enddo
          per=per+p2/ngroup
        enddo
        nu(:,:,i)=per
        per=0.d0
        do ig=1,ngroup
          do j=1,3
            call rotate_r(kap0(j,:,i),p1(j,:),u(2,ig),1)
          enddo
          do j=1,3
            call rotate_r(p1(1,j),p2(1,j),u(2,ig),1)
          enddo
          per=per+p2/ngroup
        enddo
        kap0(:,:,i)=per
      enddo
c ------------------------------------------------------------------
      allocate(seeb(3,3,0:nt))
      allocate(res(3,3,0:nt))
      do it=0,nt
        a=sig(:,:,it)
        call invers_r(3,a,3)
        call dgemm('t','n',3,3,3,1.d0,a,3,nu(1,1,it),3,0.d0,
     &             seeb(1,1,it),3)
      enddo
      seeb=seeb*evolt/2.d0  !! volt/Kelvin
      sig=sig*1.45252d12   !! Siemens/meter
      nu=nu*1.976254d13   !! [SI]
      kap0=kap0*2.6888d14   !! [SI]
      allocate(kap_el(3,3,0:nt))
      do it=0,nt
        t=tmi+dt*it
        res(:,:,it)=sig(:,:,it)
        call invers_r(3,res(1,1,it),3)
        call dgemm('n','n',3,3,3,1.d0,nu(1,1,it),3,seeb(1,1,it),3,0.d0,
     &             kap_el(1,1,it),3)
        do j=1,3
          do i=1,3
            kap_el(i,j,it)=kap0(i,j,it)-t*kap_el(i,j,it)
          enddo
        enddo
      enddo
      res=res*1.d5    !! mOm*cm
      seeb=seeb*1.d6  !! micro-volt/Kelvin
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_ChemPot')
        write(3,*)' Chemical potential in eV'
        do it=0,nt
          t=tmi+dt*it
          write(3,'(f7.1,1x,e14.7)')t,(chemp(it)-chem_pot)*evolt/2
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_Seebeck')
        write(3,*)' Seebeck coefficient in micro-Volts/Kelvin'
        do it=0,nt
          t=tmi+dt*it
          write(3,'(f7.1,6(1x,e12.5))')t,seeb(1,1,it),seeb(2,2,it),
     &                                   seeb(3,3,it),seeb(1,2,it),
     &                                   seeb(1,3,it),seeb(2,3,it)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_Conductivity')
        write(3,*)' Electric conductivity in Siemens/meter [SI]'
        do it=0,nt
          t=tmi+dt*it
          write(3,'(f7.1,6(1x,e12.5))')t,sig(1,1,it),sig(2,2,it),
     &                                   sig(3,3,it),sig(1,2,it),
     &                                   sig(1,3,it),sig(2,3,it)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_Resistivity')
        write(3,*)' Electric resistivity in mOm*cm'
        do it=0,nt
          t=tmi+dt*it
          write(3,'(f7.1,6(1x,e12.5))')t,res(1,1,it),res(2,2,it),
     &                                   res(3,3,it),res(1,2,it),
     &                                   res(1,3,it),res(2,3,it)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_El_therm_cond')
        write(3,*)' Electronic thermal conductivity in W/K/m'
        do it=0,nt
          t=tmi+dt*it
          write(3,'(f7.1,6(1x,e12.5))')t,kap_el(1,1,it),kap_el(2,2,it),
     &                                   kap_el(3,3,it),kap_el(1,2,it),
     &                                   kap_el(1,3,it),kap_el(2,3,it)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_Nu')
        do it=0,nt
          t=tmi+dt*it
          write(3,'(f7.1,6(1x,e12.5))')t,nu(1,1,it),nu(2,2,it),
     &                                   nu(3,3,it),nu(1,2,it),
     &                                   nu(1,3,it),nu(2,3,it)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_Figure_of_Merit')
        write(3,*)' Figure of merit zT'
        do it=0,nt
          t=tmi+dt*it
          s2=(seeb(1,1,it)+seeb(2,2,it)+seeb(3,3,it))/3*1.d-6
          s2=s2*s2
          s=(sig(1,1,it)+sig(2,2,it)+sig(3,3,it))/3
          kel=(kap_el(1,1,it)+kap_el(2,2,it)+kap_el(3,3,it))/3
          zt=s2*s*t/(kel+2.75d0)
          f=pi**2/3.d0*0.74258387d-8
          kel=f*s*t
          zt1=s2*s*t/(kel+2.75d0)
          write(3,'(f7.1,2(1x,e12.5))')t,zt,zt1
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_Power_Factor')
        write(3,*)' Power Factor (micro Watt/cm/K^2)'
        do it=0,nt
          t=tmi+dt*it
          s2=(seeb(1,1,it)+seeb(2,2,it)+seeb(3,3,it))/3*1.d-6
          s2=s2*s2
          s=(sig(1,1,it)+sig(2,2,it)+sig(3,3,it))/3
          f=s2*s
          write(3,'(f7.1,1x,e12.5)')t,f*1.d4
        enddo
        close(3)
      endif
      deallocate(eb,deb,sig0,sig,nu,seeb,kp,wg,ebnd,chemp,res,kap0,
     &           kap_el)
      call ending
      end
