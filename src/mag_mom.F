      SUBROUTINE mag_mom
	  use atom_mod
	  use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      IMPLICIT none
	  integer :: isort,iat0,iatom,mt,mt1,irad,ig,j,i
      real*8 :: totmom(3),spmom(3),aorbmom(3),spin(3),aorb(3),totj(3),
     &          sqpi4,sumspinx,sumspiny,sumch,g_fac,magnet(3),
     &          const,ror,charge,eff_m,s2,diskr,as,al,aj2,spinroix,
     &          spinroiy,spinroiz,chint,dqdall,ddot,al2,aj,g_e,magproj
	  real*8, allocatable :: work(:)
	  allocate(work(0:maxnrad))
      sqpi4=1.d0/sqrt(4.d0*pi)
      g_e=2.0023d0
      sumspinx=0.d0
      sumspiny=0.d0
      sumspinz=0.d0
      sumch=0.d0
	  e_b_ext=0.d0
	  magnet_o=0.d0
	  magnet_s=0.d0
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do isort=1,nsort
        iat0=iat_1(isort)
	    if(magn.eq.2) then
c     &&&&&&  The moments &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	      call orb_tot(work,totmom,spmom,aorbmom,isort)
	    endif
        do iatom=iat0,natom
          if(is(iatom).ne.isort) goto 3
	      if(maswrk)write(iun,*)' Iatom',iatom,'  ',txtel(isort)
	      spin=0.d0
	      aorb=0.d0
	      totj=0.d0
          mt=indmt(1,isort,1)
	      g_fac=2.d0
          if(irel.lt.2) then
            mt1=indmt(1,isort,nspin)
            const=dfloat(nspin)*sqpi4/2.d0
            do irad=0,nrad(isort)
              ror=ro(mt+irad)+ro(mt1+irad)
              work(irad)=ror*weight(irad,isort)*const
            enddo
            charge=dqdall(h(isort),work,nrad(isort))
            do irad=0,nrad(isort)
              ror=ro(mt1+irad)-ro(mt+irad)
              work(irad)=ror*weight(irad,isort)*const
            enddo
            spin(3)=dqdall(h(isort),work,nrad(isort))/2
            totj(3)=2*spin(3)
            eff_m=spin(3)
            as=spin(3)
            magnet_s=magnet_s-2*spin
          else if(irel.eq.2) then
            do irad=0,nrad(isort)
              ror=ro(mt+irad)
              work(irad)=ror*weight(irad,isort)*sqpi4
            enddo
            charge=dqdall(h(isort),work,nrad(isort))
	        eff_m=0.d0
            as=spin(3)
	        if(magn.eq.2) then
              do ig=1,ngroup
  		        if(ip(iatom,ig).eq.iat0) then
	              call rot_mag(ig,spmom,spin)
	              call rot_mag(ig,aorbmom,aorb)
	              call rot_mag(ig,totmom,totj)
                  magnet_s=magnet_s-spin*g_e
                  magnet_o=magnet_o-aorb
	              s2=ddot(3,spin,1,spin,1)
	              diskr=1.d0+4.d0*s2
                  as=0.5d0*(sqrt(diskr)-1.d0)
	              al2=ddot(3,aorb,1,aorb,1)
	              diskr=1.d0+4.d0*al2
                  al=0.5d0*(sqrt(diskr)-1.d0)
	              aj2=ddot(3,totj,1,totj,1)
	              diskr=1.d0+4.d0*aj2
                  aj=0.5d0*(sqrt(diskr)-1.d0)
                  if(abs(aj)<1.d-8) then
                    eff_m=0.d0
                    g_fac=1.d0
                  else
                    g_fac=1.d0+(aj*(aj+1.d0)+as*(as+1.d0)-al*(al+1.d0))
     &                        /2.d0/aj/(aj+1.d0)
	                eff_m=sqrt(aj*(aj+1.d0))*g_fac
	              endif
	              goto 2
                endif
	          enddo
2	 	      continue
            endif   !!! only for magn=2
          endif   !!! over irel
          if(maswrk) then
		    write(iun,1001)(spin(i),i=1,3),as
	        if(magn.eq.2) then
		      write(iun,1002)(aorb(i),i=1,3),al
		      write(iun,1003)(totj(i),i=1,3),aj
              write(iun,*)' Full Site Magnetization Vector:'
              magnet=-aorb-g_e*spin
              magproj=dot_product(magnet,totj)/aj2
              write(iun,1006)(magnet(i),i=1,3)
              write(iun,*)' Site Magnetization along J (m_B/u.c.):',
     &                                               magproj
	        endif
		    write(iun,1006)charge,g_fac,eff_m
            if(nspin==2.or.magn==2) then
              if(abs(spin(1))>1.d-6) print *,' For isort = ',
     &                                   isort,' M_X = ',-2*spin(1)
              if(abs(spin(2))>1.d-6) print *,' For isort = ',
     &                                   isort,' M_Y = ',-2*spin(2)
              if(abs(spin(3))>1.d-6) print *,' For isort = ',
     &                                   isort,' M_Z = ',-2*spin(3)
            endif
	      endif
          sumspinx=sumspinx+spin(1)
          sumspiny=sumspiny+spin(2)
          sumspinz=sumspinz+spin(3)
          sumch=sumch+charge
3		  continue
	    enddo  !! over iatom belonging to isort
      enddo     !!! over isort
c     &&&&&& INTERSTITIAL REGION &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      spinroix=0.d0
      spinroiy=0.d0
      spinroiz=0.d0
      chint=0.d0
      if(irel.lt.2) then
        if(.not.complex_ro) then
          do j=1,nplwro
            spinroiz=spinroiz+real(rointr(j,nspin)-rointr(j,1))*sovr(j)
            chint=chint+real(rointr(j,1)+rointr(j,nspin))*sovr(j)
          enddo
        else if(complex_ro) then
          do j=1,nplwro
            spinroiz=spinroiz+real(rointr(j,nspin)-rointr(j,1))*sovr(j)
     &          -imag(rointr(j,nspin)-rointr(j,1))*sovi(j)
            chint=chint+real(rointr(j,1)+rointr(j,nspin))*sovr(j)
     &          -imag(rointr(j,1)+rointr(j,nspin))*sovi(j)
          enddo
        endif
        chint=chint*dfloat(nspin)/2.d0*amega
        spinroiz=spinroiz*amega/2
      else if(irel.eq.2) then
        chint=ddot(nplwro,real(rointr),1,sovr,1)
        if(complex_ro) chint=chint-ddot(nplwro,imag(rointr),1,sovi,1)
        chint=chint*amega
        if(magn.eq.2) then
          if(inv.lt.2) then
            spinroix=spinroix+ddot(nplwro,real(spintr(1,:)),1,sovr,1)
            spinroiy=spinroiy+ddot(nplwro,real(spintr(2,:)),1,sovr,1)
            spinroiz=spinroiz+ddot(nplwro,real(spintr(3,:)),1,sovr,1)
          endif
          if(inv.eq.0) then
            spinroix=spinroix-ddot(nplwro,imag(spintr(1,:)),1,sovi,1)
            spinroiy=spinroiy-ddot(nplwro,imag(spintr(2,:)),1,sovi,1)
            spinroiz=spinroiz-ddot(nplwro,imag(spintr(3,:)),1,sovi,1)
          endif
        endif
        spinroix=spinroix*amega/g_e
        spinroiy=spinroiy*amega/g_e
        spinroiz=spinroiz*amega/g_e
      endif   !!! over irel
      sumspinx=sumspinx+spinroix
      sumspiny=sumspiny+spinroiy
      sumspinz=sumspinz+spinroiz
      sumch=sumch+chint
      if(nspin==2) magnet_s(3)=magnet_s(3)-2*spinroiz
      if(magn==2) then
        magnet_s(1)=magnet_s(1)-spinroix*g_e
        magnet_s(2)=magnet_s(2)-spinroiy*g_e
        magnet_s(3)=magnet_s(3)-spinroiz*g_e
      endif
      if(maswrk) then
	    write(iun,*)' Cell Sum : Spin'
        write(iun,1004)sumspinx,sumspiny,sumspinz
        if(irel==2) then
          write(iun,*)' Full Cell Magnetization Vector:'
          write(iun,1006)(magnet_s(i)+magnet_o(i),i=1,3)
        endif
        if(b_extval>1.d-8) then
          write(iun,*)' Spin susceptibility'
          as=magnet_s(3)  !sqrt(dot_product(magnet_s,magnet_s))
          as=as/amega/b_extval
          write(iun,*)' Volume: ',3.34588d-4*as
          as=2.98575d-5*as*amega
          write(iun,*)'   Mass: ',as/molar_mass,' cm^3/g'
          write(iun,*)'  Molar: ',as,' cm^3/mol'
          if(irel==2) then
            write(iun,*)' Total susceptibility'
            totmom=magnet_s+magnet_o
            as=totmom(3) !sqrt(dot_product(totmom,totmom))
            as=as/amega/b_extval
            write(iun,*)' Volume: ',3.34588d-4*as
            as=2.98575d-5*as*amega
            write(iun,*)'   Mass: ',as/molar_mass,' cm^3/g'
            write(iun,*)'  Molar: ',as,' cm^3/mol'
          endif
        endif
	  endif
1001  format(' S_x = ',f10.6,' S_y = ',f10.6,' S_z = ',f10.6,
     &       ' S = ',f10.6)
1002  format(' L_x = ',f10.6,' L_y = ',f10.6,' L_z = ',f10.6,
     &       ' L = ',f10.6)
1003  format(' J_x = ',f10.6,' J_y = ',f10.6,' J_z = ',f10.6,
     &       ' J = ',f10.6)
1004  format(' MS_x= ',f10.6,' MS_y= ',f10.6,' MS_z= ',f10.6)
1006  format(' Charge = ',f10.6,' G_factor = ',e16.6,' M_eff = ',f10.6)
	  deallocate(work)
      END
