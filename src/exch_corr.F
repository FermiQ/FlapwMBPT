      SUBROUTINE exch_corr(ROUP,RODN,ild,igga,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,VXCUP,VXCDN)
*************************************************************
*  Calculates exchange-correlation energy and potential     *
*  after LSDA or GGA91/96 [see PRL 77, 3865 (1996)].        *
*    ROUP=Rho_up                                            *
*    RODN=Rho_dn                                            *
*    GRUP=Nabla_a(Rho_up)         (  a = x,y,z or -1,0,1)   *
*    GRDN=Nabla_a(Rho_dn)         (  a = x,y,z or -1,0,1)   *
*    G2UP=Nabla_a*Nabla_b(Rho_up) (a,b = x,y,z or -1,0,1)   *
*    G2DN=Nabla_a*Nabla_b(Rho_dn) (a,b = x,y,z or -1,0,1)   *
*    ALUP=Nabla^2(Rho_up) (Laplasian)                       *
*    ALDN=Nabla^2(Rho_dn) (Laplasian)                       *
*  The ILDA=IXC-MOD(IXC,100) is avaiable with the following *
*  parametrizations:                                        *
*  ILDA=0 - no exchange-correlation                         *
*      =1 - U.von.Barth and L.Hedin,                        *
*           J.Phys.C5, 1629 (1972)                          *
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
*  To get GGA91 with specific LSDA, set IXC=ILDA+100        *
*  To get GGA96 with specific LSDA, set IXC=ILDA+200        *
*************************************************************
      use manager_mod
      IMPLICIT none
      integer, intent(in) :: ild,igga
      real*8, intent(in) :: gro(3,3),gmod(3,3),romod(3),rolap(2)
      real*8, intent(out) :: excup,excdn,vxcup,vxcdn
      integer :: i
      real*8 :: rrr,sss,roup,rodn,exldaup,vxldaup,exldadn,vxldadn,
     &          eclda,vcldaup,vcldadn,exggaup,vxggaup,exggadn,vxggadn,
     &          absg2up,grrup,grup,gmup,absgup,alup,absg2dn,grrdn,grdn,
     &          gmdn,absgdn,aldn,absg2,grdm,absg,alps,ecgga,vcggaup,
     &          vcggadn

c      SCALING>>>>>.......
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
        IF(ILD.EQ.0)RETURN
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
        EXCUP=EXCUP+EXLDAUP
        EXCDN=EXCDN+EXLDADN
        VXCUP=VXCUP+VXLDAUP
        VXCDN=VXCDN+VXLDADN
C  Local correlation:
        ECLDA=0.D0
        VCLDAUP=0.D0
        VCLDADN=0.D0
******
        IF(ILD.EQ.1.OR.ILD.EQ.2.OR.ILD.EQ.3)
     &  CALL LSDCORB(ROUP,RODN,ECLDA,VCLDAUP,VCLDADN,ILD)
        IF(ILD.EQ.4)
     &  CALL LSDCORV(ROUP,RODN,ECLDA,VCLDAUP,VCLDADN)
        IF(ILD.EQ.5)
     &  CALL LSDCORP(ROUP,RODN,ECLDA,VCLDAUP,VCLDADN)
        IF(ILD.EQ.7)
     &  CALL LSDCHL(ROUP,RODN,ECLDA,VCLDAUP,VCLDADN)
******
        EXCUP=EXCUP+ECLDA
        EXCDN=EXCDN+ECLDA
        VXCUP=VXCUP+VCLDAUP
        VXCDN=VXCDN+VCLDADN
C  *****************************************************
C  GENERALIZED GRADIENT APPROXIMATION AFTER PERDEW ET.AL
C  *****************************************************
C  Non-local exchange:
        EXGGAUP=0.D0
        VXGGAUP=0.D0
        EXGGADN=0.D0
        VXGGADN=0.D0
******
      if(igga/=0) then
        absg2up=0.d0
        grrup=0.d0
        do i=1,3
          grup=gro(i,2)
          gmup=gmod(i,2)
          absg2up=absg2up+grup*gmup
          grrup=grrup+gro(i,3)*gro(i,2)
        enddo
        absgup=romod(2)
        alup=rolap(2)
        absg2dn=0.d0
        grrdn=0.d0
        do i=1,3
          grdn=gro(i,1)
          gmdn=gmod(i,1)
          absg2dn=absg2dn+grdn*gmdn
          grrdn=grrdn+gro(i,3)*gro(i,1)
        enddo
        absgdn=romod(1)
        aldn=rolap(1)
        absg2=0.d0
        grdm=0.d0
        do i=1,3
          grdn=gro(i,3)
          gmdn=gmod(i,3)
          absg2=absg2+grdn*gmdn
          grdm=grdm+gro(i,3)*(gro(i,2)-gro(i,1))
        enddo
        absg=romod(3)
        alps=alup+aldn
        IF(IGGA==1) then
c -------------- GGA91 Exchange for spin up ------------------------
          CALL GGA91EXCH(2.D0*ROUP,2.d0*absgup,4.d0*absg2up,2.D0*ALUP,
     &                   EXGGAUP,VXGGAUP)
c ---------------- GGA91 Exchange for spin dn ----------------------
          CALL GGA91EXCH(2.D0*ROdn,2.d0*absgdn,4.d0*absg2dn,2.D0*ALdn,
     &                   EXGGAdn,VXGGAdn)
        else IF(IGGA==2) then
c -------------------- GGA96 Exchange for spin up ------------------
          CALL GGA96EXCH(2.D0*ROUP,2.d0*absgup,4.d0*absg2up,2.D0*ALUP,
     &                   EXLDAUP,VXLDAUP,EXGGAUP,VXGGAUP)
c -------------------- GGA96 Exchange for spin dn ------------------
          CALL GGA96EXCH(2.D0*RODN,2.d0*absgdn,4.d0*absg2dn,2.D0*ALDN,
     &                   EXLDADN,VXLDADN,EXGGADN,VXGGADN)
        else if(igga==3) then
          call pbe_x_sr(2.D0*roup,2.d0*absgup,4.d0*absg2up,2.D0*alup,
     &                  excup,vxcup)
          call pbe_x_sr(2.D0*rodn,2.d0*absgdn,4.d0*absg2dn,2.D0*aldn,
     &                  excdn,vxcdn)
        endif
      endif
******
      if(igga/=3) then
        EXCUP=EXCUP+EXGGAUP
        EXCDN=EXCDN+EXGGADN
        VXCUP=VXCUP+VXGGAUP
        VXCDN=VXCDN+VXGGADN
      endif
C  Non-local correlation:
      ECGGA=0.D0
      VCGGAUP=0.D0
      VCGGADN=0.D0
******
      IF(IGGA.EQ.1.AND.ILD.NE.6) !!! GGA91 Correlation
     &  CALL GGA91CORR(ROUP,RODN,absg,absg2,grdm,alps,
     &  ECLDA,VCLDAUP,VCLDADN,ECGGA,VCGGAUP,VCGGADN)
      IF(IGGA.EQ.2.AND.ILD.NE.6) !!! GGA96 Correlation
     &  CALL GGA96CORR(ROUP,RODN,absg,absg2,grrup,grrdn,alps,ECLDA,
     &                 VCLDAUP,VCLDADN,ECGGA,VCGGAUP,VCGGADN)
******
      EXCUP=EXCUP+ECGGA
      EXCDN=EXCDN+ECGGA
      VXCUP=VXCUP+VCGGAUP
      VXCDN=VXCDN+VCGGADN
      END












      SUBROUTINE GGA91EXCH(RHO,absg,absg2,ALRHO,EXGGA,VXGGA)
******************************************************************
*  Returns GGA91 correction to the exchange energy/potential.    *
*  for non-spin polarized system. For both spins call first      *
*  for rho=2*rho_up, then for rho=2*rho_dn.                      *
*  This program works with J.Perdew subroutines.                 *
******************************************************************
      use units_mod
      IMPLICIT none
      real*8 :: rho,absg,absg2,alrho,exgga,vxgga,d,fk,s,u,v,ex,vx,vx2,
     &          vxl
      EXGGA=0.D0
      VXGGA=0.D0
C Check case rho=0
      IF(RHO.LE.1.D-14)RETURN
C Prepare to call Perdew program:
      D=RHO
      FK=(3.D0*PI*PI*D)**(1.D0/3.D0)
      S=ABSG/(2.D0*FK*D)
      U=ABSG2/(D*D*(2.D0*FK)**3)
      V=ALRHO/(D*(2.D0*FK)**2)
C Call Perdew program:
        CALL PERDEXCH(D,S,U,V,EX,VX,VX2,VXL)
C Return in Rydbergs (Hartrees inside Perdew's program)
        EXGGA=EX*2.D0
        VXGGA=VX*2.D0
        END



      SUBROUTINE GGA91CORR(ROUP,RODN,absg,absg2,grdm,alps,
     &  ECLDA,VCLDAUP,VCLDADN,ECGGA,VCGGAUP,VCGGADN)
******************************************************************
*  Returns GGA91 correction to the correlation energy/potential. *
*  This program works with J.Perdew subroutines.                 *
******************************************************************
      use units_mod
      IMPLICIT none
      real*8 :: roup,rodn,absg,absg2,grdm,alps,eclda,vcldaup,vcldadn,
     &          ecgga,vcggaup,vcggadn,rho,d,zet,rs,ddrs,fk,sk,g,t,uu,vv,
     &          ww,ec,eczet,ecrs,h,vcup,vcdn,vc2up,vc2dn,vclup,vcldn
      ECGGA=0.D0
      VCGGAUP=0.D0
      VCGGADN=0.D0
C Check case rho=0
      RHO=ROUP+RODN
      IF(RHO.LE.1.D-14)RETURN
C Prepare to call Perdew program:
      D=ROUP+RODN
      ZET=(ROUP-RODN)/D
      RS=(3.D0/(4.D0*PI*D))**(1.D0/3.D0)
      DDRS=-9.D0/(4.D0*PI*RS**4)
      FK=(3.D0*PI*PI*D)**(1.D0/3.D0)
      SK=SQRT(4.D0*FK/PI)
      G=((1.D0+ZET)**(2.D0/3.D0)+(1.D0-ZET)**(2.D0/3.D0))/2.D0
      T=ABSG/(2.D0*G*SK*D)
      UU=ABSG2/(D*D*(2.D0*SK*G)**3)
      VV=ALPS/(D*(2.D0*SK*G)**2)
      WW=(GRDM/D-ZET*ABSG**2/D)/(D*(2*SK*G)**2)
C EC=Correlation energy in LDA, ECZET=dEC/dZET(rho=const),
C ECRS=dEC/dRS(zet=const.!); Hartree units now:
        EC=ECLDA*0.5D0
        ECZET=0.25D0*(VCLDAUP-ECLDA)-0.25D0*(VCLDADN-ECLDA)
        ECRS=0.25D0*((VCLDAUP-ECLDA)*(1.D0+ZET)+
     &               (VCLDADN-ECLDA)*(1.D0-ZET))/D*DDRS
C Call Perdew program:
        CALL PERDCORR(RS,ZET,T,UU,VV,WW,H,VCUP,VCDN,
     &  VC2UP,VC2DN,VCLUP,VCLDN,FK,SK,G,EC,ECRS,ECZET)
C Return in Rydbergs (Hartrees inside Perdew's program)
        ECGGA=H*2.D0
        VCGGAUP=VCUP*2.D0
        VCGGADN=VCDN*2.D0
        END





      SUBROUTINE PERDEXCH(D,S,U,V,EX,VX,VX2,VXL)
******************************************************************
*  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM       *
*  INPUT D : DENSITY                                             *
*  INPUT S:  ABS(GRAD D)/(2*KF*D)                                *
*  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)       *
*  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)                          *
*  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX) *
*  OUTPUT:  CONTRIBUTIONS TO VX FROM U AND V (added)             *
******************************************************************
      IMPLICIT none
      real*8 :: d,s,u,v,ex,vx,vx2,vxl,a1,a2,a3,a4,ax,a,b1,thrd,thrd4,
     &          fac,s2,s3,s4,p0,p1,p2,p3,p4,f,p5,p6,p7,p8,p9,fs,fss,p10,
     &          p11
        DATA A1,A2,A3,A4/0.19645D0,0.27430D0,0.15084D0,100.D0/
        DATA AX,A,B1/-0.738558766382022D0,7.7956D0,0.004D0/
        DATA THRD,THRD4/0.33333333333333D0,1.33333333333333D0/
        FAC = AX*D**THRD
        S2 = S*S
        S3 = S2*S
        S4 = S2*S2
        P0 = 1.D0/SQRT(1.D0+A*A*S2)
        P1 = LOG(A*S+1.D0/P0)
        P2 = EXP(-A4*S2)
        P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
        P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
        F = P3*P4
        EX = FAC*F-FAC ! without local contribution
C  LOCAL EXCHANGE OPTION
C       EX = FAC
C  ENERGY DONE. NOW THE POTENTIAL:
        P5 = B1*S2-(A2-A3*P2)
        P6 = A1*S*(P1+A*S*P0)
        P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
        FS = P3*(P3*P5*P6+P7)
        P8 = 2.D0*S*(B1-A3*A4*P2)
        P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
        P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
        P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
        FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
        VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)-FAC*THRD4 !without local part
C  CONTRIBUTIONS TO VX FROM SECOND-ORDER GRADIENTS (ADDED)
        VX2 = -FAC*FSS
        VXL = -FAC*FS
C  LOCAL EXCHANGE OPTION:
C       VX = FAC*THRD4
        RETURN
        END




      SUBROUTINE PERDCORR(RS,ZET,T,UU,VV,WW,H,DVCUP,DVCDN,
     &  VC2UP,VC2DN,VCLUP,VCLDN,FK,SK,G,EC,ECRS,ECZET)
******************************************************************
*  GGA91 CORRELATION                                             *
*  INPUT RS: SEITZ RADIUS                                        *
*  INPUT ZET: RELATIVE SPIN POLARIZATION                         *
*  INPUT T: ABS(GRAD D)/(D*2.*KS*G)                              *
*  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)     *
*  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)                     *
*  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2               *
*  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON    *
*  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS *
*  OUTPUT VC2UP,VC2DN:  CONTRIBUTIONS TO VC FROM UU (added)      *
*  OUTPUT VCLUP,VCLDN:  CONTRIBUTIONS TO VC FROM VV (added)      *
******************************************************************
      IMPLICIT none
      real*8 :: rs,zet,t,uu,vv,ww,h,dvcup,dvcdn,vc2up,vc2dn,vclup,vcldn,
     &          fk,sk,g,ec,ecrs,eczet,xnu,cc0,cx,alf,c1,c2,c3,c4,c5,c6,
     &          a4,thrdm,thrd2,bet,delt,g3,g4,pon,b,b2,t2,t4,t6,rs2,rs3,
     &          q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3,h0,h1,r4,gz,fact5,
     &          fac,bg,bec,q8,q9,q82,h0b,h0rs,fact0,fact1,h0bt,h0rst,
     &          h0z,h0t,h0zt,fact2,fact3,h0tt,h1rs,fact4,h1rst,h1z,h1t,
     &          ccrs,rsthrd,h1zt,h1tt,hrs,hrst,ht,htt,hz,hzt,comm,pref
      DATA XNU,CC0,CX,ALF/15.75592D0,0.004235D0,-0.001667212D0,0.09D0/
      DATA C1,C2,C3,C4/0.002568D0,0.023266D0,7.389D-6,8.723D0/
      DATA C5,C6,A4/0.472D0,7.389D-2,100.D0/
      DATA THRDM,THRD2/-0.33333333333333D0,0.66666666666667D0/
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(EXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
C  LOCAL CORRELATION OPTION:
C       H = 0.0D0
C  ENERGY DONE. NOW THE POTENTIAL:
      CCRS = (C2+2.D0*C3*RS)/Q7 -
     &  Q6*(C4+2.D0*C5*RS+3.D0*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      Q82 = Q8*Q8
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*(Q4*Q5*FACT0-DELT*FACT1)/Q82
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.D0*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-Q9*FACT3/Q82)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
C  CONTRIBUTIONS TO VC FROM SECOND-ORDER GRADIENTS (ADDED)
      VC2UP = -HTT
      VC2DN = -HTT
      VCLUP = -HT
      VCLDN = -HT
C  LOCAL CORRELATION OPTION:
C       DVCUP = 0.0D0
C       DVCDN = 0.0D0
      RETURN
      END


      SUBROUTINE GGA96EXCH(RHO,absg,absg2,ALRHO,EXLDA,VXLDA,
     &                       EXGGA,VXGGA)
******************************************************************
*  Returns GGA96 correction to the exchange energy/potential.    *
*  for non-spin polarized system. For both spins call first      *
*  for rho=2*rho_up, then for rho=2*rho_dn.                      *
******************************************************************
      use units_mod
      IMPLICIT none
      real*8 :: rho,absg,absg2,alrho,exlda,vxlda,exgga,vxgga,ckap,cmuu,
     &          fk,ss,dss,gss,dgss,fxs1,dfxs,gfxs,dgfxs,g2fxs,vsurf
C  GGA constants
        DATA CKAP/0.804000000000000000D+00/  !
     &      ,CMUU/0.219516451220895831D+00/  !  betta*pi^2/3
        EXGGA=0.D0
        VXGGA=0.D0
C  Check case rho=0
        IF(RHO.LE.1.D-14)RETURN
        IF(absg.LE.1.D-14)RETURN
C  Set electron gas data
        FK=(3.D0*PI*PI*RHO)**(1.D0/3.D0)
        SS=ABSG/(2.D0*FK*RHO)
C  Set derivatives
        DSS=-2.D0*ABSG/(3.D0*FK*RHO*RHO)
        GSS=1.D0/(2.D0*FK*RHO)
        DGSS=-2.D0/(3.D0*FK*RHO*RHO)
C  Compute GGA96 to exchange
        FXS1=CKAP-CKAP/(1.D0+CMUU*SS*SS/CKAP)
        DFXS=2.D0*CMUU*SS/(1.D0+CMUU*SS*SS/CKAP)**2*DSS
        GFXS=2.D0*CMUU*SS/(1.D0+CMUU*SS*SS/CKAP)**2*GSS
        DGFXS=2.D0*CMUU*SS/(1.D0+CMUU*SS*SS/CKAP)**2*DGSS+
     &  2.D0*CMUU*DSS/(1.D0+CMUU*SS*SS/CKAP)**2*GSS-8.D0*
     &  CMUU**2*SS**2*DSS/CKAP/(1.D0+CMUU*SS*SS/CKAP)**3*GSS
        G2FXS=
     &  2.D0*CMUU*GSS/(1.D0+CMUU*SS*SS/CKAP)**2*GSS-8.D0*
     &  CMUU**2*SS**2*GSS/CKAP/(1.D0+CMUU*SS*SS/CKAP)**3*GSS
        EXGGA=EXLDA*FXS1
        VSURF=VXLDA*GFXS*ABSG+EXLDA*DGFXS*RHO*ABSG+
     &  EXLDA*G2FXS*RHO*ABSG2/ABSG+EXLDA*GFXS*
     &  RHO/ABSG*(ALRHO-ABSG2/ABSG)
        VXGGA=VXLDA*FXS1+RHO*EXLDA*DFXS-VSURF
        END



      SUBROUTINE GGA96CORR(ROUP,RODN,absg,absg2,grrup,grrdn,alps,
     &  ECLDA,VCLDAUP,VCLDADN,ECGGA,VCGGAUP,VCGGADN)
******************************************************************
*  Returns GGA96 correction to the correlation energy/potential. *
******************************************************************
      use units_mod
      IMPLICIT none
      real*8, intent(in) :: roup,rodn,absg,absg2,grrup,grrdn,alps,eclda,
     &                      vcldaup,vcldadn
      real*8, intent(out) :: ecgga,vcggaup,vcggadn
      integer :: ispin
      real*8 :: a0,e2,cgam,cbet,rho,yy,fk,sk,ph,tt,dfk,dsk,dph,dtt,gtt,
     &          dgtt,cd,dyy,declda,dcd,ca,dca,ct1,dct1,gct1,dgct1,g2ct1,
     &          ct2,dct2,gct2,dgct2,g2ct2,ct,dct,gct,dgct,g2ct,cl,dcl,
     &          gcl,dgcl,g2cl,decgga,gecgga,dgecgga,g2ecgga
C  Fundamential constants
      DATA A0/1.D0/ ! Bohr radius
     &    ,E2/2.D0/ ! El.charge^2
C  GGA constants
      DATA CGAM/0.310906908696549009E-01/  ! (1-ln(2))/pi^2
     &    ,CBET/0.066725000000000000D+00/  !
      ECGGA=0.D0
      VCGGAUP=0.D0
      VCGGADN=0.D0
C  Check case rho=0
      RHO=ROUP+RODN
      IF(RHO.LE.1.D-14)RETURN
      IF(absg.LE.1.D-14)RETURN
C  Set electron gas data
      YY=(ROUP-RODN)/RHO
      FK=(3.D0*PI*PI*RHO)**(1.D0/3.D0)
      SK=SQRT(4.D0*FK/PI/A0)
      PH=((1.D0+YY)**(2.D0/3.D0)+(1.D0-YY)**(2.D0/3.D0))/2.D0
      TT=ABSG/(2.D0*PH*SK*RHO)
      VCGGAUP=0.D0
      VCGGADN=0.D0
      DO ISPIN=1,2                                   !!! over spins
C  Set derivatives
        IF(ISPIN.EQ.1)DYY=-(ROUP-RODN)/RHO/RHO+1.D0/RHO
        IF(ISPIN.EQ.2)DYY=-(ROUP-RODN)/RHO/RHO-1.D0/RHO
        DFK=(3.D0*PI*PI)**(1.D0/3.D0)*1.D0/3.D0*RHO**(-2.D0/3.D0)
        DSK=DFK/SQRT(PI*FK*A0)
        DPH=((1.D0+YY)**(-1.D0/3.D0)-(1.D0-YY)**(-1.D0/3.D0))/
     &  3.D0*DYY
        DTT=-ABSG/(2.D0*PH*PH*SK*SK*RHO*RHO)*
     &  (SK*RHO*DPH+SK*PH+PH*RHO*DSK)
        GTT=1.D0/(2.D0*PH*SK*RHO)
        DGTT=-1.D0/(2.D0*PH*PH*SK*SK*RHO*RHO)*
     &  (SK*RHO*DPH+SK*PH+PH*RHO*DSK)
C  Compute GGA96 to correlation
C  1)
        CD=-1.D0/(CGAM*E2/A0)*ECLDA/PH**3
        IF(ISPIN.EQ.1)DECLDA=(VCLDAUP-ECLDA)/RHO
        IF(ISPIN.EQ.2)DECLDA=(VCLDADN-ECLDA)/RHO
        DCD=-1.D0/(CGAM*E2/A0)*(DECLDA*PH-ECLDA*3.D0*DPH)/PH**4
C  2)
        CA=CBET/CGAM/(EXP(CD)-1.D0)
        DCA=-CBET/CGAM/(EXP(CD)-1.D0)**2*EXP(CD)*DCD
C  3)
        CT1=1.D0+CA*TT**2
        DCT1=2.D0*TT*DTT*CA+DCA*TT*TT
        GCT1=2.D0*TT*GTT*CA
        DGCT1=2.D0*DTT*GTT*CA+2.D0*TT*DGTT*CA+2.D0*TT*GTT*DCA
        G2CT1=2.D0*GTT*GTT*CA
        CT2=CT1+CA**2*TT**4
        DCT2=DCT1+4.D0*TT**3*DTT*CA**2+2.D0*CA*DCA*TT**4
        GCT2=GCT1+4.D0*TT**3*GTT*CA**2
        DGCT2=DGCT1+4.D0*3.D0*TT**2*DTT*GTT*CA**2+
     &  4.D0*TT**3*DGTT*CA**2+4.D0*TT**3*GTT*2.D0*CA*DCA
        G2CT2=G2CT1+4.D0*3.D0*TT**2*GTT*GTT*CA**2
        CT=CT1/CT2
        DCT=(CT2*DCT1-DCT2*CT1)/(CT2*CT2)
        GCT=(CT2*GCT1-GCT2*CT1)/(CT2*CT2)
        DGCT=(DCT2*GCT1+CT2*DGCT1-DGCT2*CT1-GCT2*DCT1)/(CT2*CT2)-
     &  2.D0*(CT2*GCT1-GCT2*CT1)/CT2**3*DCT2
        G2CT=(GCT2*GCT1+CT2*G2CT1-G2CT2*CT1-GCT2*GCT1)/(CT2*CT2)-
     &  2.D0*(CT2*GCT1-GCT2*CT1)/CT2**3*GCT2
C  4)
        CL=1.D0+CBET/CGAM*TT*TT*CT
        DCL=CBET/CGAM*(2.D0*TT*DTT*CT+TT*TT*DCT)
        GCL=CBET/CGAM*(2.D0*TT*GTT*CT+TT*TT*GCT)
        DGCL=CBET/CGAM*(2.D0*DTT*GTT*CT+2.D0*TT*DGTT*CT+
     &  2.D0*TT*GTT*DCT+2.D0*TT*DTT*GCT+TT*TT*DGCT)
        G2CL=CBET/CGAM*(2.D0*GTT*GTT*CT+2.D0*TT*GTT*GCT+
     &  2.D0*TT*GTT*GCT+TT*TT*G2CT)
C  Answer
        ECGGA=E2/A0*CGAM*PH**3*LOG(CL)
        DECGGA=E2/A0*CGAM*(3.D0*PH**2*DPH*LOG(CL)+PH**3/CL*DCL)
        GECGGA=E2/A0*CGAM*PH**3/CL*GCL
        DGECGGA=E2/A0*CGAM*(3.D0*PH**2*DPH/CL*GCL-
     &  PH**3/CL**2*DCL*GCL+PH**3/CL*DGCL)
        G2ECGGA=E2/A0*CGAM*PH**3*(G2CL/CL-GCL*GCL/CL**2)
        IF(ISPIN.EQ.1)
     &  VCGGAUP=VCGGAUP+ECGGA+RHO*DECGGA-
     &  ABSG*GECGGA-RHO*ABSG2/ABSG*G2ECGGA-
     &  RHO/ABSG*(ALPS-ABSG2/ABSG)*GECGGA
        IF(ISPIN.EQ.2)
     &  VCGGADN=VCGGADN+ECGGA+RHO*DECGGA-
     &  ABSG*GECGGA-RHO*ABSG2/ABSG*G2ECGGA-
     &  RHO/ABSG*(ALPS-ABSG2/ABSG)*GECGGA
        IF(ISPIN.EQ.1)
     &  VCGGAUP=VCGGAUP-RHO/ABSG*GRRUP*DGECGGA
        IF(ISPIN.EQ.1)
     &  VCGGADN=VCGGADN-RHO/ABSG*GRRUP*DGECGGA
        IF(ISPIN.EQ.2)
     &  VCGGAUP=VCGGAUP-RHO/ABSG*GRRDN*DGECGGA
        IF(ISPIN.EQ.2)
     &  VCGGADN=VCGGADN-RHO/ABSG*GRRDN*DGECGGA
      ENDDO                                          !!! over spins
      END
