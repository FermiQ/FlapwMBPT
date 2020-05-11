      SUBROUTINE EXCVXC(ROUP,RODN,ild0,EXCUP,EXCDN,VXCUP,VXCDN)
*************************************************************
*  Calculates exchange-correlation energy and potential     *
*  after LSDA or GGA91/96 [see PRL 77, 3865 (1996)].        *
*    ROUP=Rho_up                                            *
*    RODN=Rho_dn                                            *
*  The ILDA=IXC-MOD(IXC,100) is avaiable with the following *
*  parametrizations:                                        *
*  ILD0=0 - no exchange-correlation                         *
*      =1 - U.von.Barth and L.Hedin,                        *
*           J.Phys.C 5, 1629 (1972)                         *
*      =2 - O.E.Gunnarsson and S.Lundqvist                  *
*           Phys.Rev.B                                      *
*      =3 - V.L.Moruzzi, J.F.Janak, and A.R.Williams,       *
*           Calculated Electronic Properties of Metals      *
*           (New York, Pergamon Press, 1978)                *
*      =4 - S.H.Vosko, L.Wilk, and M.Nusair                 *
*           Can.J.Phys.58, 1200 (1980)                      *
*           parametrization of the Ceperly-Alder form:      *
*           D.M.Ceperly, B.J.Alder,                         *
*           Phys.Rev.Lett.45, 566 (1980)                    *
*      =5   J.Perdew, X.Wang (1991) programs of LSDA=4      *
*      =6 - Original Gaspar-Kohn-Sham (rho**1/3)            *
*     In case ILD0<0 only correlation is taken              *
*************************************************************
      use etot_mod
      use manager_mod
      IMPLICIT none
	  integer, intent(in) :: ild0
      real*8, intent(out) :: excup,excdn,vxcup,vxcdn
	  integer :: ild
	  real*8 :: rrr,sss,roup,rodn,exldaup,vxldaup,exldadn,vxldadn,
     &          eclda,vcldaup,vcldadn

c	   SCALING>>>>>.......
      ild=iabs(ild0)
	  rrr=roup+rodn
	  sss=scal_spin*(roup-rodn)
	  roup=0.5d0*(rrr+sss)
	  rodn=0.5d0*(rrr-sss)
c>>>>>>>>>>>>>>>>>>....................

C  Energies and potentials
      EXCUP=0.D0
      EXCDN=0.D0
      VXCUP=0.D0
      VXCDN=0.D0
      IF(ild.EQ.0)RETURN
C  ********************************
C  LOCAL SPIN DENSITY APPROXIMATION
C  ********************************
C  Local exchange
        EXLDAUP=0.D0
        VXLDAUP=0.D0
        EXLDADN=0.D0
        VXLDADN=0.D0
******
        CALL LSDEXCH(2.D0*ROUP,EXLDAUP,VXLDAUP)
        CALL LSDEXCH(2.D0*RODN,EXLDADN,VXLDADN)
******
        if(ild0>=0) then
          EXCUP=EXCUP+EXLDAUP
          EXCDN=EXCDN+EXLDADN
          VXCUP=VXCUP+VXLDAUP
          VXCDN=VXCDN+VXLDADN
        endif
C  Local correlation:
        ECLDA=0.D0
        VCLDAUP=0.D0
        VCLDADN=0.D0
******
        IF(ild.EQ.1.OR.ild.EQ.2.OR.ild.EQ.3)
     &  CALL LSDCORB(ROUP,RODN,ECLDA,VCLDAUP,VCLDADN,ild)
        IF(ild.EQ.4)
     &  CALL LSDCORV(ROUP,RODN,ECLDA,VCLDAUP,VCLDADN)
        IF(ild.EQ.5)
     &  CALL LSDCORP(ROUP,RODN,ECLDA,VCLDAUP,VCLDADN)
        IF(ild.EQ.7)
     &  CALL LSDCHL(ROUP,RODN,ECLDA,VCLDAUP,VCLDADN)
******
        EXCUP=EXCUP+ECLDA
        EXCDN=EXCDN+ECLDA
        VXCUP=VXCUP+VCLDAUP
        VXCDN=VXCDN+VCLDADN
      END


      SUBROUTINE LSDEXCH(RHO,EX,VX)
********************************************************
*  Calculates Exchange Energy and Potential            *
*  after Gaspar-Kohn-Sham (rho**1/3) for non-spin      *
*  polarized system. To take into account both spins:  *
*  Ex=rho_up*ex(2*rho_up)+rho_dn*ex(2*rho_dn).         *
*  where <<ex>> is returned by this subroutine.        *
*  Correspondingly, for the potential:                 *
*  Vx_up=vx(2*rho_up), Vx_dn=vx(2*rho_dn)              *
********************************************************
	  use units_mod
      IMPLICIT none
	  real*8 :: alf,rho,ex,vx,rs
      DATA ALF/0.521061761197848D0/ !!! (4/9/pi)^(1/3)
C Check for zero
      IF(RHO.LE.1.D-44)THEN
        EX=0.D0
        VX=0.D0
      ELSE
        RS=(3.D0/(4.D0*PI*RHO))**(1.D0/3.D0)
        EX=-3.D0/(2.D0*PI*ALF*RS)
        VX=-2.D0/(     PI*ALF*RS)
      ENDIF
      END


      SUBROUTINE LSDCORB(ROUP,RODN,EC,VCUP,VCDN,I)
********************************************************
*  Calculates Correlation Energy and potential after   *
*  Barth-Hedin.                                        *
*    I=1 - U.von.Barth and L.Hedin,                    *
*          J.Phys.C5, 1629 (1972)                      *
*    I=2 - O.E.Gunnarsson and S.Lundqvist              *
*          Phys.Rev.B                                  *
*    I=3 - V.L.Moruzzi, J.F.Janak, and A.R.Williams,   *
*          Calculated Electronic Properties of Metals  *
*          (New York, Pergamon Press, 1978)            *
********************************************************
	  use units_mod
      IMPLICIT none
	  integer :: i
      real*8 :: CCP(3),CCF(3),RP(3),RF(3),ecx,x,decx,f,y,g1,roup,rodn,
     &	      rho,ec,vcup,vcdn,rs,xp,xf,g2,dec,o
C Instrinsic formulae:
      ECX(X)=-((1.D0+X**3)*LOG(1.D0+1.D0/X)+X/2.D0-X*X-1.D0/3.D0)
      DECX(X)=-(3.D0*X*X*LOG(1.D0+1.D0/X)-(1.D0+X**3)/
     & (X*X+X)+1.D0/2.D0-2.D0*X)
        F(Y)=((1.D0+Y)**(4.D0/3.D0)+G1-2.D0)/(2.D0*(
     &  2.D0**(1.D0/3.D0)-1.D0))
C Correlation parameters:
      DATA CCP,CCF/0.0504D0,0.0666D0,0.045D0,0.0254D0
     &      ,0.0406D0,0.0357D0/RP,RF/30.D0,11.4D0,21.D0
     &      ,75.D0,15.9D0,26.5D0/
      IF(ROUP.LT.0D0)ROUP=1.D-10
      IF(RODN.LT.0D0)RODN=1.D-10
      RHO=ROUP+RODN
C Check for zero
      IF(RHO.LE.1.D-44)THEN
        EC=0.D0
        VCUP=0.D0
        VCDN=0.D0
      ELSE
        RS=(3.D0/(4.D0*PI*(ROUP+RODN)))**(1.D0/3.D0)
        XP=RS/RP(I)
        XF=RS/RF(I)
        Y=(ROUP-RODN)/(ROUP+RODN)
        G1=(1.D0-Y)**(4.D0/3.D0)
        G2=(1.D0-Y)**(1.D0/3.D0)
        XP=RS/RP(I)
        XF=RS/RF(I)
C Correlation energy
        EC=CCP(I)*ECX(XP)+(CCF(I)*ECX(XF)-CCP(I)*ECX(XP))*F(Y)
        DEC=(CCF(I)*ECX(XF)-CCP(I)*ECX(XP))*((1.D0+Y)**
     &  (1.D0/3.D0)-G2)*(4.D0/(6.D0*(2.D0**(1.D0/3.D0)-1.D0)))
C Correlation potential
        O=-1.D0
        VCUP=EC-(RS/3.D0)*(CCP(I)/RP(I)*DECX(XP)+
     & (CCF(I)/RF(I)*DECX(XF)-CCP(I)/RP(I)*
     &  DECX(XP))*F(Y))-O*(1.D0+O*Y)*DEC
        O=+1.D0
        VCDN=EC-(RS/3.D0)*(CCP(I)/RP(I)*DECX(XP)+
     & (CCF(I)/RF(I)*DECX(XF)-CCP(I)/RP(I)*
     &  DECX(XP))*F(Y))-O*(1.D0+O*Y)*DEC
      ENDIF
      END


      SUBROUTINE LSDCORV(ROUP,RODN,EC,VCUP,VCDN)
********************************************************
*  Calculates Exchange-Correlation Energy              *
*  and potential after:                                *
*        - S.H.Vosko, L.Wilk, and M.Nusair             *
*          Can.J.Phys.58, 1200 (1980)                  *
*          parametrization of the Ceperly-Alder form:  *
*          D.M.Ceperly, B.J.Alder,                     *
*          Phys.Rev.Lett.45, 566 (1980)                *
********************************************************
	  use units_mod
      IMPLICIT none
      real*8 :: oth,fth,aa,ap,xp0,bp,cp,qp,cp1,cp2,cp3,af,xf0,bf,cf,qf,
     &          cf1,cf2,cf3,roup,rodn,rho,ec,vcup,vcdn,rho1,rho2,rs1,rs,
     &          x,xpx,s,sp,sm,s4,fs,beta,atnp,ecp,xfx,atnf,ecf,dfs,
     &          dbeta,tp1,tf1,ucp,ucf,uc0,uc20,uc10,duc,duc1,duc2,uc1,
     &          uc2
      DATA OTH/0.33333333333333D0/,FTH/1.33333333333333D0/
      DATA AA/0.519842D0/
      DATA AP,XP0,BP,CP,QP,CP1,CP2,CP3/0.0621814D0,-0.10498D0
     &      ,3.72744D0,12.9352D0,6.1519908D0,1.2117833D0,1.1435257D0
     &      ,-0.031167608D0/
      DATA AF,XF0,BF,CF,QF,CF1,CF2,CF3/0.0310907D0,-0.32500D0,
     &                                 7.060428D0,18.0578D0,4.7309269D0,
     &                                 2.9847935D0,2.7100059D0,
     &                                 -0.1446006D0/
C      **********************************
      IF(ROUP.LT.0D0)ROUP=1.D-10
      IF(RODN.LT.0D0)RODN=1.D-10
      RHO=ROUP+RODN
C Check for zero
      IF(RHO.LE.1.D-44)THEN
        EC=0.D0
        VCUP=0.D0
        VCDN=0.D0
      ELSE
        RHO1=ROUP
        RHO2=RODN
        RS1=(4.D0*PI*RHO/3.D0)**OTH
        RS=1.D0/RS1
C Correlation Energy
        X=SQRT(RS)
        XPX=X*X+BP*X+CP
        S=(RHO2-RHO1)/RHO
        SP=1.D0+S
        SM=1.D0-S
        S4=S**4-1.D0
        FS=(SP**FTH+SM**FTH-2.D0)/AA
        BETA=1.D0/(2.74208D0+3.182D0*X+0.09873D0*X*X+0.18268D0*X**3)
        ATNP=ATAN(QP/(2.D0*X+BP))
        ECP=AP*(LOG(X*X/XPX)+CP1*ATNP-CP3*(LOG((X-XP0)**2/XPX)
     &  +CP2*ATNP))
        XFX=X*X+BF*X+CF
        ATNF=ATAN(QF/(2.D0*X+BF))
        ECF=AF*(LOG(X*X/XFX)+CF1*ATNF-CF3*(LOG((X-XF0)**2/XFX)
     &  +CF2*ATNF))
        EC=ECP+FS*(ECF-ECP)*(1.D0+S4*BETA)
C Correlation Potential
        DFS=FTH*(SP**OTH-SM**OTH)/AA
        DBETA=-(0.27402*X+0.09873+1.591/X)*BETA**2
        TP1=(X*X+BP*X)/XPX
        TF1=(X*X+BF*X)/XFX
        UCP=ECP-AP/3.D0*(1.D0-TP1-CP3*(X/(X-XP0)-TP1-XP0*X/XPX))
        UCF=ECF-AF/3.D0*(1.D0-TF1-CF3*(X/(X-XF0)-TF1-XF0*X/XFX))
        UC0=UCP+(UCF-UCP)*FS
        UC20=UC0+(ECF-ECP)*SM*DFS
        UC10=UC0-(ECF-ECP)*SP*DFS
        DUC=(UCF-UCP)*BETA*S4*FS+(ECF-ECP)*(-RS/3.)*DBETA*S4*FS
        DUC2=DUC+(ECF-ECP)*BETA*SM*(4.D0*S**3*FS+S4*DFS)
        DUC1=DUC-(ECF-ECP)*BETA*SP*(4.D0*S**3*FS+S4*DFS)
        UC1=UC10+DUC1
        UC2=UC20+DUC2
        VCUP=UC1
        VCDN=UC2
      ENDIF
      END




      SUBROUTINE LSDCORP(ROUP,RODN,EC,VCUP,VCDN)
*****************************************************
*  Uniform-gas correlation of Perdew and Wang 1991  *
*  Practically identical to Vosko et.al.            *
*****************************************************
	  use units_mod
      IMPLICIT none
      real*8 :: gam,fzz,thrd,thrd4,roup,rodn,ec,vcup,vcdn,rho,rs,zet,f,
     &          z4,eu,eurs,ep,eprs,alfm,alfrsm,ecrs,fz,eczet,comm
      DATA GAM,FZZ/0.5198421D0,1.709921D0/
      DATA THRD,THRD4/0.33333333333333D0,1.33333333333333D0/
C  INPUT
      IF(ROUP.LT.0D0)ROUP=1.D-10
      IF(RODN.LT.0D0)RODN=1.D-10
      RHO=ROUP+RODN
C Check for zero
      IF(RHO.LE.1.D-44)THEN
        EC=0.D0
        VCUP=0.D0
        VCDN=0.D0
        RETURN
      ENDIF
      RS=(3.D0/(4.D0*PI*(ROUP+RODN)))**(1.D0/3.D0)
      ZET=(ROUP-RODN)/(ROUP+RODN)
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      CALL GCORP(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     &    0.49294D0,1.00D0,RS,EU,EURS)
      CALL GCORP(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     &    0.62517D0,1.00D0,RS,EP,EPRS)
      CALL GCORP(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     &    0.49671D0,1.00D0,RS,ALFM,ALFRSM)
C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
c       ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     &        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
C  RETURN IN RYDBERG'S
      EC=2.D0*EC
      VCUP=2.D0*VCUP
      VCDN=2.D0*VCDN
      END
      SUBROUTINE GCORP(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
C  CALLED BY SUBROUTINE LSDCORP
      IMPLICIT none
      real*8 :: a,a1,b1,b2,b3,b4,p,rs,gg,ggrs,p1,q0,rs12,rs32,rsp,q1,q2,
     &          q3
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      END
      subroutine LSDCHL(roup,rodn,ec,vcup,vcdn)
********************************************************
*  calculates exchange-correlation potential           *
*     =7 - L.Hedin and S.Lundqvist                     *
********************************************************
	  use units_mod
      implicit none
	  real*8 :: roup,rodn,ec,vcup,vcdn,alf,rs,cp,rpp,zp
      parameter(alf=0.521062d0)
      if(roup.lt.0d0)roup=1.d-15
      if(rodn.lt.0d0)rodn=1.d-15
      rs=(3.d0/(4.d0*pi*(roup+rodn)))**(1.d0/3.d0)
      CP=0.038608D0
      RPp=24.3D0
      ZP=RS/RPp
      VCup=-CP*DLOG(1.0D0+1.0D0/ZP)
      VCdn=vcup
      ec=-Cp*((1.0D0+Zp**3)*DLOG(1.0D0+1.0D0/Zp)+Zp/2.0D0
     *            -Zp*Zp-1.d0/3.d0)
      end
