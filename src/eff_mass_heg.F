      subroutine eff_mass_heg_1
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
      integer :: i,i0,n,k_0
      real*8 :: pi2a,z0,s0,ds0,meff,zk,hh,w,dedk,fm,f0,fp,aaa,bbb,ccc,
     &          k_f,meff1,f1,x0,x1,x2,s2,f2
      real*8, allocatable :: f(:),a(:),b(:),c(:),d(:)
      if(n_kl_fermi-1<1) return
      if(n_kl_fermi+1>n_line) return
      pi2a=(pi+pi)/par
      allocate(f(n_line))
      allocate(a(n_line-1))
      allocate(b(n_line-1))
      allocate(c(n_line-1))
      allocate(d(n_line-1))
      hh=(qabs_line(2)-qabs_line(1))*pi2a
      if(ubi=='dft') then
        k_f=sqrt(chem_pot-vxc_heg(nspin))
      else
c -------- Spline approximation of Sigma(w=0) -------------------
      do i=1,n_line
        i0=kline_in_npnt(i)
        f(i)=sigx_k_heg(i0,nspin)
        if(ubi/=' hf') f(i)=f(i)+sig_0_heg(i0,nspin)
      enddo
      call spline_hmg(hh,f,a,b,c,d,n_line-1,0,0.d0,0.d0)
c ---------- Effective Fermi momentum -----------------------------
        k_0=0
        f0=f(1)-chem_pot
        do i=2,n_line
          f1=(hh*(i-1))**2+f(i)-chem_pot
          if(f0*f1<0.d0) then
            k_0=i-1
            goto 1
          endif
          f0=f1
        enddo
1       x0=0.d0
        x1=hh
2       x2=0.5d0*(x0+x1)
        s2=a(k_0)+b(k_0)*x2+c(k_0)*x2**2+d(k_0)*x2**3
        f2=(hh*(k_0-1)+x2)**2+s2-chem_pot
        if(f1*f2>0.d0) then
          f1=f2
          x1=x2
        else
          f0=f2
          x0=x2
        endif
        if(x1-x0>1.d-9) goto 2
        x2=0.5d0*(x0+x1)
        k_f=hh*(k_0-1)+x2
      endif
      n=k_f/hh
      w=k_f-hh*n
      if(ubi=='dft') then
        s0=vxc_heg(nspin)
        ds0=0.d0
        z0=1.d0
      else
c -------- Sig_0 at k_fermi ------------------------------------------
        s0=a(k_0)+b(k_0)*x2+c(k_0)*x2**2+d(k_0)*x2**3
        ds0=b(k_0)+2.d0*c(k_0)*x2+3.d0*d(k_0)*x2**2
c -------- Z at k_fermi ---------------------------------------------
      do i=1,n_line
        i0=kline_in_npnt(i)
        f(i)=z_ren_heg(i0,nspin)
        if(ubi==' hf') f(i)=1.d0
      enddo
      call spline_hmg(hh,f,a,b,c,d,n_line-1,0,0.d0,0.d0)
        z0=a(k_0)+b(k_0)*x2+c(k_0)*x2**2+d(k_0)*x2**3
      endif
c --------- Effective mass m*/m --------------------------------------
      zk=1.d0+ds0/2/k_f
      meff=z0/zk
c ----Effective mass from dE/dK -----------------------------------
      fm=spectrum(n_kl_fermi-1)
      f0=spectrum(n_kl_fermi)
      fp=spectrum(n_kl_fermi+1)
      ccc=f0
      bbb=(fp-fm)/2/hh
      aaa=(fp-bbb*hh-ccc)/hh/hh
      dedk=2.d0*aaa*w+bbb
      meff1=k_f/dedk*2.d0
c -------------------------------------------------------------------      
      if(maswrk) then
        write(iun,*)' K_Fermi = ',k_fermi
        write(iun,*)' K_Fermi_effective = ',k_f
        write(iun,*)' Z at K_Fermi = ',1.d0/z0
        write(iun,*)' ZK at K_Fermi = ',zk
        write(iun,*)' Effective mass ratio M*/M = ',meff
        write(iun,*)' M*/M from de/dk = ',meff1
      endif
      deallocate(f,a,b,c,d)
      end
      
      
      subroutine eff_mass_heg1
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
      integer :: i,k_0
      real*8 :: pi2a,z0,s0,ds0,meff,zk,w,dedk,f0,k_f,meff1,f1,x0,x1,x2,
     &          s2,f2
      real*8, allocatable :: f(:),a(:),b(:),c(:),d(:)
      pi2a=(pi+pi)/par
      allocate(f(0:n_k_heg))
      allocate(a(n_k_heg))
      allocate(b(n_k_heg))
      allocate(c(n_k_heg))
      allocate(d(n_k_heg))
      if(ubi=='dft') then
        k_f=sqrt(chem_pot-vxc_heg(nspin))
        w=x_fermi+asinh((k_f-k_fermi)/alf_k_heg)
        k_0=w/step_k_x
        x2=0.d0
      else
c -------- Spline approximation of Sigma(w=0) -------------------
        do i=0,n_k_heg
          f(i)=sigx_k_heg(i,nspin)
          if(ubi/=' hf') f(i)=f(i)+sig_0_heg(i,nspin)
        enddo
        call spline_inhmg(k_m_heg,f,a,b,c,d,n_k_heg,0,0.d0,0.d0)
c ---------- Effective Fermi momentum -----------------------------
        k_0=0
        f0=f(0)-chem_pot
        do i=1,n_k_heg
          f1=k_m_heg(i)**2+f(i)-chem_pot
          if(f0*f1<0.d0) then
            k_0=i
            goto 1
          endif
          f0=f1
        enddo
1       x0=0.d0
        x1=k_m_heg(k_0)-k_m_heg(k_0-1)
2       x2=0.5d0*(x0+x1)
        s2=a(k_0)+b(k_0)*x2+c(k_0)*x2**2+d(k_0)*x2**3
        f2=(k_m_heg(k_0-1)+x2)**2+s2-chem_pot
        if(f1*f2>0.d0) then
          f1=f2
          x1=x2
        else
          f0=f2
          x0=x2
        endif
        if(x1-x0>1.d-9) goto 2
        x2=0.5d0*(x0+x1)
        k_f=k_m_heg(k_0-1)+x2
      endif
      w=k_f-k_m_heg(k_0-1)
      if(ubi=='dft') then
        s0=vxc_heg(nspin)
        ds0=0.d0
        z0=1.d0
      else
c -------- Sig_0 at k_fermi ------------------------------------------
        s0=a(k_0)+b(k_0)*x2+c(k_0)*x2**2+d(k_0)*x2**3
        ds0=b(k_0)+2.d0*c(k_0)*x2+3.d0*d(k_0)*x2**2
c -------- Z at k_fermi ---------------------------------------------
        do i=0,n_k_heg
          f(i)=z_ren_heg(i,nspin)
          if(ubi==' hf') f(i)=1.d0
        enddo
        call spline_inhmg(k_m_heg,f,a,b,c,d,n_k_heg,0,0.d0,0.d0)
        z0=a(k_0)+b(k_0)*x2+c(k_0)*x2**2+d(k_0)*x2**3
      endif
c --------- Effective mass m*/m --------------------------------------
      zk=1.d0+ds0/2/k_f
      meff=z0/zk
c ----Effective mass from dE/dK -----------------------------------
      do i=0,n_k_heg
        f(i)=spectrum(i)
      enddo
      call spline_inhmg(k_m_heg,f,a,b,c,d,n_k_fermi+10,0,0.d0,0.d0)
      dedk=b(k_0)+2.d0*c(k_0)*x2+3.d0*d(k_0)*x2**2
      meff1=k_f/dedk*2.d0
c -------------------------------------------------------------------      
      if(maswrk) then
        write(iun,*)' K_Fermi = ',k_fermi
        write(iun,*)' K_Fermi_effective = ',k_f
        write(iun,*)' Z at K_Fermi = ',1.d0/z0
        write(iun,*)' ZK at K_Fermi = ',zk
        write(iun,*)' Effective mass ratio M*/M = ',meff
        write(iun,*)' M*/M from de/dk = ',meff1
      endif
      deallocate(f,a,b,c,d)
      end