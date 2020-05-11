	  subroutine heg
	  use atom_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      nloc=iexch/100
      dim_heg=3
      ro_el_gas=3.d0/4.d0/pi/rs_el_gas**3
      amega=1.d0/ro_el_gas
      par=amega**(1.d0/3.d0)
      k_fermi=(3.d0*pi**2*ro_el_gas)**(1.d0/3.d0)
      emindos=-max(2*k_fermi**2,k_fermi**2+1.5d0)
      emaxdos=max(2*k_fermi**2,k_fermi**2+1.d0)
      nu_plasmon=sqrt(16.d0*pi*ro_el_gas)
      if(maswrk) then
        write(iun,*)' K_Fermi = ',k_fermi
        write(iun,*)' Electron Gas Density = ',ro_el_gas
        write(iun,*)' E_kin = ',0.6d0*k_fermi**2
        write(iun,*)' E_x = ',-1.5d0*k_fermi/pi
        write(iun,*)' Plasmon frequency = ',nu_plasmon
        write(iun,*)' Non-interacting compressibility = ',
     &              k_fermi/2.d0/pi**2/ro_el_gas**2
c        call compress_from_etot
c        call compress_from_etot_gw
      endif
      omega_geom=omega_geom*k_fermi**2*6
      omega_max=omega_max*k_fermi**2*6
      nu_geom=nu_geom*k_fermi**2*6
      nu_max=nu_max*k_fermi**2*6
      call tau_frequency_meshes_spl
      if(dim_heg==2) allocate(u_opt(2,2,48))
      if(dim_heg==3) allocate(u_opt(3,3,48))
      call group_heg
      numvol=50
      if(iter_psi+iter_bsp/=0) call heg_mesh_c
      call heg_rad
      end
      
      
      
      subroutine compress_from_etot
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer :: i,n
      real*8 :: hh,et(0:20),rs(0:20),f(0:20),mu(0:20),p00(0:20),
     &          dens(0:20),spl(20,4),kp1,kp0(0:20),kf(0:20),
     &          vrtq(0:20),vrtnu(0:20),zf(0:20),mst(0:20),kpe(0:20),
     &          kpm(0:20),kpp(0:20),kpv(0:20),f0s(0:20)
      if(maswrk) then
        open(3,file='c:\RESEARCH\QFT\RESULT\HEG\Plots\F_QMC')
        read(3,*)
        n=12
        do i=0,n
          read(3,*)rs(i),et(i),mu(i),dens(i),kf(i),kp0(i)
        enddo
        close(3)
c ------- Compressibility from Etot ------------------------------
        kpe=0.d0
        do i=0,n
          f(i)=dens(i)*et(i)
        enddo
        call spline_inhmg(dens,f,spl(1,1),spl(1,2),spl(1,3),spl(1,4),n,
     &                    0,0.d0,0.d0)
        do i=1,n-1
          kp1=dens(i)**2*2*spl(i+1,3)
          kpe(i)=1.d0/kp1
        enddo
c ------- Compressibility from Chemical Potential -----------------
        call spline_inhmg(rs,mu,spl(1,1),spl(1,2),spl(1,3),spl(1,4),n,
     &                    0,0.d0,0.d0)
        do i=0,n-1
          kp1=-spl(i+1,2)/4.d0/pi/rs(i)**2
          kpm(i)=1.d0/kp1
        enddo
        hh=rs(n)-rs(n-1)
        kp1=-(spl(n,2)+2*spl(n,3)*hh+3*spl(n,4)*hh*hh)/4.d0/pi/rs(n)**2
        kpm(n)=1.d0/kp1
        open(3,
     &   file='c:\RESEARCH\QFT\RESULT\HEG\Plots\compressibility.aqmc')
        do i=0,n
          write(3,'(f5.2,4(1x,e12.5))')rs(i),kpe(i)/kp0(i),
     &                                       kpm(i)/kp0(i)
        enddo
        close(3)
        open(3,file='c:\RESEARCH\QFT\RESULT\HEG\Plots\F_D')
        read(3,*)
        n=11
        do i=0,n
          read(3,*)rs(i),et(i),mu(i),dens(i),kf(i),kp0(i),p00(i),
     &             vrtnu(i),vrtq(i),zf(i),mst(i)
        enddo
        close(3)
c ------- Compressibility from Etot ------------------------------
        kpe=0.d0
        do i=0,n
          f(i)=dens(i)*et(i)
        enddo
        call spline_inhmg(dens,f,spl(1,1),spl(1,2),spl(1,3),spl(1,4),n,
     &                    0,0.d0,0.d0)
        do i=1,n-1
          kp1=dens(i)**2*2*spl(i+1,3)
          kpe(i)=1.d0/kp1
        enddo
c ------- Compressibility from Chemical Potential -----------------
        call spline_inhmg(rs,mu,spl(1,1),spl(1,2),spl(1,3),spl(1,4),n,
     &                    0,0.d0,0.d0)
        do i=0,n-1
          kp1=-spl(i+1,2)/4.d0/pi/rs(i)**2
          kpm(i)=1.d0/kp1
        enddo
        hh=rs(n)-rs(n-1)
        kp1=-(spl(n,2)+2*spl(n,3)*hh+3*spl(n,4)*hh*hh)/4.d0/pi/rs(n)**2
        kpm(n)=1.d0/kp1
c -------- Compressibility from the Polarizability -----------------
        do i=0,n
          kpp(i)=-p00(i)/dens(i)**2
        enddo
c -------- Compressibility from the Vertices -----------------
        do i=0,n
          kpv(i)=kp0(i)*vrtq(i)/vrtnu(i)*mst(i)
          f0s(i)=vrtnu(i)/vrtq(i)-1.d0
        enddo
        open(3,
     &   file='c:\RESEARCH\QFT\RESULT\HEG\Plots\compressibility.d')
        do i=0,n
          write(3,'(f5.2,5(1x,e12.5))')rs(i),kpe(i)/kp0(i),
     &                                       kpm(i)/kp0(i),
     &                                       kpp(i)/kp0(i),
     &                                       kpv(i)/kp0(i),f0s(i)
        enddo
        close(3)
      endif
      end
      
      
      
      subroutine compress_from_etot_gw
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer :: i,n
      real*8 :: hh,et(0:90),rs(0:90),f(0:90),mu(0:90),
     &          dens(0:90),spl(90,4),kp1,kp0(0:90),kpe(0:90),
     &          kpm(0:90)
      real*8, allocatable :: v(:)
      if(maswrk) then
        open(3,file='c:\RESEARCH\QFT\RESULT\HEG\Plots\F_A0')
        read(3,*)
        n=17
        do i=0,n
          read(3,*)rs(i),et(i),mu(i),dens(i),kp0(i)
        enddo
        close(3)
c ------- Compressibility from Etot ------------------------------
        allocate(v(0:n))
        kpe=0.d0
        do i=0,n
          f(i)=dens(i)*et(i)
          v(i)=4.d0*pi*rs(i)**3/3.d0
        enddo
        call spline_inhmg(dens,f,spl(1,1),spl(1,2),spl(1,3),spl(1,4),n,
     &                    0,0.d0,0.d0)
        do i=1,n-1
          kp1=dens(i)**2*2*spl(i+1,3)
          kpe(i)=1.d0/kp1
        enddo
        call spline_inhmg(v,et,spl(1,1),spl(1,2),spl(1,3),spl(1,4),n,
     &                    0,0.d0,0.d0)
        do i=1,n-1
          kp1=v(i)*2*spl(i+1,3)
c          kpe(i)=1.d0/kp1
        enddo
        deallocate(v)
c ------- Compressibility from Chemical Potential -----------------
        call spline_inhmg(rs,mu,spl(1,1),spl(1,2),spl(1,3),spl(1,4),n,
     &                    0,0.d0,0.d0)
        do i=0,n-1
          kp1=-spl(i+1,2)/4.d0/pi/rs(i)**2
          kpm(i)=1.d0/kp1
        enddo
        hh=rs(n)-rs(n-1)
        kp1=-(spl(n,2)+2*spl(n,3)*hh+3*spl(n,4)*hh*hh)/4.d0/pi/rs(n)**2
        kpm(n)=1.d0/kp1
        open(3,
     &   file='c:\RESEARCH\QFT\RESULT\HEG\Plots\compressibility.a0')
        do i=0,n
          write(3,'(f5.2,5(1x,e12.5))')rs(i),kpe(i)/kp0(i),
     &                                       kpm(i)/kp0(i)
        enddo
        close(3)
      endif
      end
