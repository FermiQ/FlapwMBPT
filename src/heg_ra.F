      subroutine heg_ra(rs,ksi)
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(in) :: rs,ksi
	integer :: iq,i_nu,k
	real*8 :: q2,t0,ex,ec,dec_rs,dec_ksi,d2ec_rs,d2ec_ksi,
     &	        vxc(2),tc,lams0,lamn0,lama0,lamsinf,g0,gams,laman,w,as,
     &          bs,cs,g01,alf1,alf,q8,q6,pi8,lamninf,lamainf,gama,gamn,
     &          om_f,pi2a,aa,ca,bet,ba,an,cn,bn,q4
      pi2a=(pi+pi)/par
      om_f=(k_fermi)**2
      pi8=8.d0*pi
	call heg_qmc_prm(rs,ksi,t0,ex,ec,dec_rs,dec_ksi,d2ec_rs,d2ec_ksi,
     &	               vxc,tc)
      laman=1.d0-1.5d0*(2.d0*pi/3.d0)**(2.d0/3.d0)*rs*d2ec_ksi
      lamn0=-0.11d0*rs/(1.d0+0.33d0*rs)*laman
      alf=(4.d0/9.d0/pi)**(1.d0/3.d0)
      lams0=1.d0+pi/3.d0*alf*rs*rs*dec_rs-pi/6.d0*alf*rs**3*d2ec_rs
     &     -lamn0
      lama0=laman-lamn0
      lamsinf=0.6d0-2.d0*pi*alf/5.d0*(rs*rs*dec_rs+2.d0*rs*ec)
      lamninf=3.d0*pi*alf*rs*(ec+rs*dec_rs)
      lamainf=0.d0
      g0=0.5d0*(1.d0+2*0.193d0*rs)
     &  /(1.d0+0.525d0*rs+0.193d0*0.525d0*rs*rs)**2
      alf1=0.9d0
      gams=9.d0/16.d0*(1.d0-g0)*lamsinf+0.25d0+0.75d0*(1.d0-1.d0/alf1)
      gama=9.d0/8.d0*lamainf+0.25d0
      gamn=0.68d0
      g01=1.d0-g0
      if(maswrk) then
        write(iun,*)' Results from QMC parametrization :'
        write(iun,*)' T0 = ',t0
        write(iun,*)' TC = ',tc
        write(iun,*)' EX = ',ex
        write(iun,*)' EC = ',ec
        write(iun,*)' dec_drs = ',dec_rs
        write(iun,*)' dec_dksi = ',dec_ksi
        write(iun,*)' d2ec_drs2 = ',d2ec_rs
        write(iun,*)' d2ec_dksi2 = ',d2ec_ksi
        write(iun,*)' Lambda_s_inf = ',lamsinf
        write(iun,*)' g(0) = ',g0
      endif
	do i_nu=0,n_nu
	  w=w_nu(i_nu)/4.d0/om_f
	  as=lamsinf+(lams0-lamsinf)/(1.d0+gams**2*w**2)
	  cs=0.75d0*lamsinf/g01
	  cs=cs-(4.d0/3.d0-1.d0/alf1+cs)/(1.d0+gams*w)
	  bs=as/(3*as*(1.d0+w)**4-8.d0*g01*(1.d0+w)**3/3.d0
     &	                       -2.d0*cs*g01*(1.d0+w)**4)
c        write(iun,'(4(1x,e15.8))')w,as,bs,cs
c ------ Constant and quadratic terms are here -------------------  
	  do iq=1,nstar
	    q2=(pi2a*qabs_star(iq)/2/k_fermi)**2
	    q6=q2**3
	    q8=q2**4
	    g_fac_s(iq,i_nu)=(as*q2+2.d0*bs*g01*q8/3.d0)
     &	                    /(1.d0+cs*q2+bs*q8)
	  enddo
	  aa=lamainf+(lama0-lamainf)/(1.d0+gama**2*w**2)
	  ca=1.5d0*lamainf-(1.d0/3.d0+1.5d0*lamainf)/(1.d0+gama**2*w)
	  bet=1.d0/3.d0*(4.d0*g0-1.d0)
     & 	   +1.d0/3.d0*(1.d0-2.d0*g0)*gama**2*w*w/(1.d0+gama**2*w*w)
	  ba=aa/(3*aa*(1.d0+w)**4-4.d0*bet*(1.d0+w)**3
     &	                       -3.d0*ca*bet*(1.d0+w)**4)
c        write(iun,'(4(1x,e15.8))')w,aa,ba,ca
c ------ Constant and quadratic terms are here -------------------  
	  do iq=1,nstar
	    q2=(pi2a*qabs_star(iq)/2/k_fermi)**2
	    q6=q2**3
	    q8=q2**4
	    g_fac_a(iq,i_nu)=-1.d0/3.d0*(1.d0-2.d0*g0)
     &	    *gama**2*w*w/(1.d0+gama**2*w*w)
     &      +(aa*q2+ba*bet*q8)/(1.d0+ca*q2+ba*q8)
	  enddo
	  an=lamninf+(lamn0-lamninf)/(1.d0+gamn**2*w**2)
	  cn=3.d0*gamn*w/1.18d0/(1.d0+gamn*w)
	  cn=cn-((lamn0+lamninf/3.d0)/(lamn0+2*lamninf/3.d0)+
     &	  3.d0*gamn*w/1.18d0/(1.d0+gamn*w))/(1.d0+gamn**2*w*w)
	  bn=-1.5d0/lamninf/(1.d0+gamn*w)**2
     &	  *(an+lamninf+2.d0*lamninf*cn*(1.d0+gamn*w)/3.d0+
     &    sqrt((an+lamninf+2.d0*lamninf*cn*(1.d0+gamn*w)/3.d0)**2
     &        +4.d0*an*lamninf/3.d0))
c        write(iun,'(4(1x,e15.8))')w,an,bn,cn
c ------ Constant and quadratic terms are here -------------------  
	  do iq=1,nstar
	    q2=(pi2a*qabs_star(iq)/2/k_fermi)**2
	    q4=q2*q2
	    q6=q2**3
	    q8=q2**4
	    g_fac_n(iq,i_nu)=(an*q2-bn*lamninf/3.d0*q6)
     &	                    /(1.d0+cn*q2+bn*q4)
	  enddo
	enddo
      do k=1,nstar
        do i_nu=0,n_nu
          f_xc(k,i_nu)=-v_q_heg(k)*(g_fac_s(k,i_nu)+g_fac_n(k,i_nu))
        enddo
      enddo
      f_xc(1,:)=f_xc(2,:)
      end
      
      
      SUBROUTINE heg_qmc_prm(rs,ksi,t0,ex,ec,dec_rs,dec_ksi,d2ec_rs,
     &                       d2ec_ksi,vxc,tc)
	use units_mod
      IMPLICIT none
      real*8, intent(in) :: rs,ksi
      real*8, intent(out) :: t0,ex,ec,dec_rs,dec_ksi,d2ec_rs,d2ec_ksi,
     &                       vxc(2),tc
	real*8 :: gam,fzz,thrd,thrd4,f,z4,eu,eurs,ep,eprs,alfm,alfrsm,
     &	        ecrs,fz,eczet,comm,eu2rs,ep2rs,al2rs,d2z
c      DATA GAM,FZZ/0.5198421D0,1.709921D0/
      thrd=1.d0/3.d0
      thrd4=4.d0/3.d0
      gam=2.d0**thrd4-2.d0
      fzz=8.d0/9.d0/gam
      F = ((1.D0+ksi)**THRD4+(1.D0-ksi)**THRD4-2.D0)/GAM
      CALL Gparam(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     &            0.49294D0,1.00D0,RS,EU,EURS,eu2rs)
      CALL Gparam(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     &            0.62517D0,1.00D0,RS,EP,EPRS,ep2rs)
      CALL Gparam(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     &            0.49671D0,1.00D0,RS,ALFM,ALFRSM,al2rs)
C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
c       ALFC = -ALFM
      Z4 = ksi**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ksi)**THRD-(1.D0-ksi)**THRD)/GAM
      ECZET = 4.D0*(ksi**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     &        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ksi*ECZET
      vxc(1) = COMM + ECZET
      vxc(2) = COMM - ECZET
c -------- Kinetic energy ----------------------------------
      t0=0.15d0/rs/rs*(9.d0*pi/4.d0)**(2.d0/3.d0)
     &  *((1.d0+ksi)**(5.d0/3.d0)+(1.d0-ksi)**(5.d0/3.d0))
      tc=-4.d0*ec+1.5d0*((1.d0+ksi)*vxc(1)+(1.d0-ksi)*vxc(2))
c -------- Exchange energy ----------------------------------
      ex=-1.5d0/4.d0/pi/rs*(9.d0*pi/4.d0)**thrd
     &  *((1.d0+ksi)**thrd4+(1.d0-ksi)**thrd4)
c -------- First derivatives ----------------------------
      dec_rs=ecrs
      dec_ksi=eczet
c -------- Second derivatives ----------------------------
      d2ec_rs=EU2RS*(1.D0-F*Z4)+EP2RS*F*Z4-AL2rs*F*(1.D0-Z4)/FZZ
      d2z = 4.d0/9.d0*((1.D0+ksi)**(THRD-1.d0)
     &                +(1.D0-ksi)**(THRD-1.d0))/GAM
      d2ec_ksi=4.D0*(3.d0*ksi**2*F+ksi**3*fz)*(EP-EU+ALFM/FZZ)
     &        +d2z*(Z4*(EP-EU)-(1.D0-Z4)*ALFM/FZZ)
     &        +fz*(4.d0*ksi**3*(EP-EU)+4.d0*ksi**3*ALFM/FZZ)
C  RETURN IN RYDBERG'S
      EC=2.D0*EC
      vxc=2.d0*vxc
      t0=2.d0*t0
      tc=2.d0*tc
      ex=2.d0*ex
      dec_rs=2.d0*dec_rs
      dec_ksi=2.d0*dec_ksi
      d2ec_rs=2.d0*d2ec_rs
      d2ec_ksi=2.d0*d2ec_ksi
      end
      SUBROUTINE Gparam(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS,gg2rs)
      IMPLICIT none
      real*8, intent(in) :: A,A1,B1,B2,B3,B4,P,RS
      real*8, intent(out) :: GG,GGRS,gg2rs
	real*8 :: p1,q0,rs12,rs32,rsp,q1,q2,q3,q0d,q1dd,t
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = rs*RS12
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      q0d=-2.d0*a*a1
      q1dd = A*(-0.5d0*B1/RS32+1.5D0*B3/RS12+2.D0*B4*p*P1*RS**(P-1.d0))
      t=(q0d*q3+q0*q1dd)*(q1*q1+q1)
      t=t-q0*q3*q3*(2.d0*q1+1.d0)
      t=t/q1/q1/(1.d0+q1)**2
      gg2rs=-q0d*q3/q1/(1.d0+q1)-t
      END
      
      
      
      subroutine heg1_ra(rs,ksi)
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(in) :: rs,ksi
	integer :: iq,i_nu
	real*8 :: q2,t0,ex,ec,dec_rs,dec_ksi,d2ec_rs,d2ec_ksi,
     &	        vxc(2),tc,lams0,lamn0,lama0,lamsinf,g0,gams,laman,w,as,
     &          bs,cs,g01,alf1,alf,q8,q6,pi8,lamninf,lamainf,gama,gamn,
     &          om_f,aa,ca,bet,ba,an,cn,bn,q4
      om_f=(k_fermi)**2
      pi8=8.d0*pi
	call heg_qmc_prm(rs,ksi,t0,ex,ec,dec_rs,dec_ksi,d2ec_rs,d2ec_ksi,
     &	               vxc,tc)
      laman=1.d0-1.5d0*(2.d0*pi/3.d0)**(2.d0/3.d0)*rs*d2ec_ksi
      lamn0=-0.11d0*rs/(1.d0+0.33d0*rs)*laman
      alf=(4.d0/9.d0/pi)**(1.d0/3.d0)
      lams0=1.d0+pi/3.d0*alf*rs*rs*dec_rs-pi/6.d0*alf*rs**3*d2ec_rs
     &     -lamn0
      lama0=laman-lamn0
      lamsinf=0.6d0-2.d0*pi*alf/5.d0*(rs*rs*dec_rs+2.d0*rs*ec)
      lamninf=3.d0*pi*alf*rs*(ec+rs*dec_rs)
      lamainf=0.d0
      g0=0.5d0*(1.d0+2*0.193d0*rs)
     &  /(1.d0+0.525d0*rs+0.193d0*0.525d0*rs*rs)**2
      alf1=0.9d0
      gams=9.d0/16.d0*(1.d0-g0)*lamsinf+0.25d0+0.75d0*(1.d0-1.d0/alf1)
      gama=9.d0/8.d0*lamainf+0.25d0
      gamn=0.68d0
      g01=1.d0-g0
      if(maswrk) then
        write(iun,*)' Results from QMC parametrization :'
        write(iun,*)' T0 = ',t0
        write(iun,*)' TC = ',tc
        write(iun,*)' EX = ',ex
        write(iun,*)' EC = ',ec
        write(iun,*)' dec_drs = ',dec_rs
        write(iun,*)' dec_dksi = ',dec_ksi
        write(iun,*)' d2ec_drs2 = ',d2ec_rs
        write(iun,*)' d2ec_dksi2 = ',d2ec_ksi
        write(iun,*)' Lambda_s_inf = ',lamsinf
        write(iun,*)' g(0) = ',g0
      endif
	do i_nu=0,n_nu
	  w=w_nu(i_nu)/4.d0/om_f
	  as=lamsinf+(lams0-lamsinf)/(1.d0+gams**2*w**2)
	  cs=0.75d0*lamsinf/g01
	  cs=cs-(4.d0/3.d0-1.d0/alf1+cs)/(1.d0+gams*w)
	  bs=as/(3*as*(1.d0+w)**4-8.d0*g01*(1.d0+w)**3/3.d0
     &	                       -2.d0*cs*g01*(1.d0+w)**4)
c        write(iun,'(4(1x,e15.8))')w,as,bs,cs
c ------ Constant and quadratic terms are here -------------------  
	  do iq=0,n_q_heg
	    q2=(q_m_heg(iq)/2/k_fermi)**2
	    q6=q2**3
	    q8=q2**4
	    g_fac_s(iq,i_nu)=(as*q2+2.d0*bs*g01*q8/3.d0)
     &	                    /(1.d0+cs*q2+bs*q8)
	  enddo
	  aa=lamainf+(lama0-lamainf)/(1.d0+gama**2*w**2)
	  ca=1.5d0*lamainf-(1.d0/3.d0+1.5d0*lamainf)/(1.d0+gama**2*w)
	  bet=1.d0/3.d0*(4.d0*g0-1.d0)
     & 	   +1.d0/3.d0*(1.d0-2.d0*g0)*gama**2*w*w/(1.d0+gama**2*w*w)
	  ba=aa/(3*aa*(1.d0+w)**4-4.d0*bet*(1.d0+w)**3
     &	                       -3.d0*ca*bet*(1.d0+w)**4)
c        write(iun,'(4(1x,e15.8))')w,aa,ba,ca
c ------ Constant and quadratic terms are here -------------------  
	  do iq=0,n_q_heg
	    q2=(q_m_heg(iq)/2/k_fermi)**2 
	    q6=q2**3
	    q8=q2**4
	    g_fac_a(iq,i_nu)=-1.d0/3.d0*(1.d0-2.d0*g0)
     &	    *gama**2*w*w/(1.d0+gama**2*w*w)
     &      +(aa*q2+ba*bet*q8)/(1.d0+ca*q2+ba*q8)
	  enddo
	  an=lamninf+(lamn0-lamninf)/(1.d0+gamn**2*w**2)
	  cn=3.d0*gamn*w/1.18d0/(1.d0+gamn*w)
	  cn=cn-((lamn0+lamninf/3.d0)/(lamn0+2*lamninf/3.d0)+
     &	  3.d0*gamn*w/1.18d0/(1.d0+gamn*w))/(1.d0+gamn**2*w*w)
	  bn=-1.5d0/lamninf/(1.d0+gamn*w)**2
     &	  *(an+lamninf+2.d0*lamninf*cn*(1.d0+gamn*w)/3.d0+
     &    sqrt((an+lamninf+2.d0*lamninf*cn*(1.d0+gamn*w)/3.d0)**2
     &        +4.d0*an*lamninf/3.d0))
c        write(iun,'(4(1x,e15.8))')w,an,bn,cn
c ------ Constant and quadratic terms are here -------------------  
	  do iq=0,n_q_heg
	    q2=(q_m_heg(iq)/2/k_fermi)**2 
	    q4=q2*q2
	    q6=q2**3
	    q8=q2**4
	    g_fac_n(iq,i_nu)=(an*q2-bn*lamninf/3.d0*q6)
     &	                    /(1.d0+cn*q2+bn*q4)
	  enddo
	enddo
      do iq=1,n_q_heg
        do i_nu=0,n_nu
          f_xc(iq,i_nu)=-v_q_heg(iq)*(g_fac_s(iq,i_nu)+g_fac_n(iq,i_nu))
        enddo
      enddo
      f_xc(0,:)=f_xc(1,:)
	call output_heg1_f_xc
      call output_heg1_g_fac
      end