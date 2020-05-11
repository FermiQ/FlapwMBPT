      subroutine inp_gga(sintet,costet,sinphi,cosphi,irad,yl,d12yl,
     &                   rotot,romod,gro,gmod,rolap,drodr,d2rodr2,dspdr,
     &                   d2spdr2,isort,totm,sm)
      use atom_mod
      use manager_mod
      use solid_mod
      implicit none
c     d12yl(1) = d Y / d teta
c     d12yl(2) = d Y / d fi
c     d12yl(3) = d2 Y / d teta2
c     d12yl(4) = d2 Y / d teta / d fi
c     d12yl(5) = d2 Y / d fi2
      integer, intent(in) :: irad,isort
      real*8, intent(in) :: sintet,costet,sinphi,cosphi,yl(limlpb),
     &                      d12yl(limlpb,5),
     &                      drodr(0:nrad(isort),nsym(isort),nspin),
     &                      d2rodr2(0:nrad(isort),nsym(isort),nspin),
     &                      dspdr(0:nrad(isort),*),
     &                      d2spdr2(0:nrad(isort),*)
      real*8, intent(out) :: rotot(2),gro(3,3),gmod(3,3),romod(3),
     &                       rolap(2)
      integer :: ndimv,ispin,lm,lget,isym,mt,i,isymb,mtb
      real*8 :: g(3),rorr(2),ror(2),rot(2),rof(2),rort(2),
     &          rorf(2),rott(2),rotf(2),roff(2),dmdr(3),
     &          dmdt(3),dmdf(3),d2mdr2(3),d2mdt2(3),d2mdf2(3),
     &          d2mdrdt(3),d2mdrdf(3),d2mdtdf(3),totm(3),
     &          r1,r2,rst,rst2,rv,fl,ro0,ro1,ro2,sm,ddot,dsmdr,
     &          dsmdt,dsmdf,d2smdr2,d2smdrdt,d2smdrdf,d2smdt2,d2smdf2,
     &          d2smdtdf,sm_lap,rrr,rrt,rrf,rtt,rtf,rff,rvv,rt,rf
      rv=r(irad,isort)
      r1=1.d0/rv
      r2=r1**2
      rst=r1/sintet
      rst2=rst**2
      ndimv=(lmpb(isort)+1)**2
      if(nspin.eq.2) then
        do ispin=1,nspin
          ror(ispin)=0.d0
          rot(ispin)=0.d0
          rof(ispin)=0.d0
          rorr(ispin)=0.d0
          rort(ispin)=0.d0
          rorf(ispin)=0.d0
          rott(ispin)=0.d0
          rotf(ispin)=0.d0
          roff(ispin)=0.d0
          rolap(ispin)=0.d0
          rotot(ispin)=0.d0
          do lm=1,ndimv
            if(sym(lm,isort)) then
              fl=dfloat(lget(lm))
              fl=fl*(fl+1.d0)
              isym=lmsym(lm,isort)
              mt=indmt(isym,isort,ispin)
              ro0=ro(mt+irad)
              ro1=drodr(irad,isym,ispin)
              ro2=d2rodr2(irad,isym,ispin)
              rotot(ispin)=rotot(ispin)+ro0*yl(lm)
              rolap(ispin)=rolap(ispin)
     &                  +(ro2+2.d0*r1*ro1-fl*r2*ro0)*yl(lm)
              ror(ispin)=ror(ispin)+ro1*yl(lm)
              rot(ispin)=rot(ispin)+ro0*d12yl(lm,1)
              rof(ispin)=rof(ispin)+ro0*d12yl(lm,2)
              rorr(ispin)=rorr(ispin)+ro2*yl(lm)
              rort(ispin)=rort(ispin)+ro1*d12yl(lm,1)
              rorf(ispin)=rorf(ispin)+ro1*d12yl(lm,2)
              rott(ispin)=rott(ispin)+ro0*d12yl(lm,3)
              rotf(ispin)=rotf(ispin)+ro0*d12yl(lm,4)
              roff(ispin)=roff(ispin)+ro0*d12yl(lm,5)
            endif
          enddo  !!! over lm
        enddo
      else if(nspin.eq.1) then
        ror(1)=0.d0
        rot(1)=0.d0
        rof(1)=0.d0
        rorr(1)=0.d0
        rort(1)=0.d0
        rorf(1)=0.d0
        rott(1)=0.d0
        rotf(1)=0.d0
        roff(1)=0.d0
        rolap(1)=0.d0
        rotot(1)=0.d0
        do lm=1,ndimv
          if(sym(lm,isort)) then
            fl=dfloat(lget(lm))
            fl=fl*(fl+1.d0)
            isym=lmsym(lm,isort)
            mt=indmt(isym,isort,1)
            ro0=ro(mt+irad)/2.d0
            ro1=drodr(irad,isym,1)/2.d0
            ro2=d2rodr2(irad,isym,1)/2.d0
            rotot(1)=rotot(1)+ro0*yl(lm)
            rolap(1)=rolap(1)+(ro2+2.d0*r1*ro1-fl*r2*ro0)*yl(lm)
            ror(1)=ror(1)+ro1*yl(lm)
            rot(1)=rot(1)+ro0*d12yl(lm,1)
            rof(1)=rof(1)+ro0*d12yl(lm,2)
            rorr(1)=rorr(1)+ro2*yl(lm)
            rort(1)=rort(1)+ro1*d12yl(lm,1)
            rorf(1)=rorf(1)+ro1*d12yl(lm,2)
            rott(1)=rott(1)+ro0*d12yl(lm,3)
            rotf(1)=rotf(1)+ro0*d12yl(lm,4)
            roff(1)=roff(1)+ro0*d12yl(lm,5)
          endif
        enddo  !!! over lm
        if(magn.eq.1) then
c ------------ Copying to ispin=2 ------------------------------
          rotot(2)=rotot(1)
          rolap(2)=rolap(1)
          ror(2)=ror(1)
          rot(2)=rot(1)
          rof(2)=rof(1)
          rorr(2)=rorr(1)
          rort(2)=rort(1)
          rorf(2)=rorf(1)
          rott(2)=rott(1)
          rotf(2)=rotf(1)
          roff(2)=roff(1)
c -----------------------------------------------------------
        else if(magn.eq.2) then
c ----------- Calculate derivatives of Magnetization -------------------
          do i=1,3
            dmdr(i)=0.d0
            dmdt(i)=0.d0
            dmdf(i)=0.d0
            d2mdr2(i)=0.d0
            d2mdrdt(i)=0.d0
            d2mdrdf(i)=0.d0
            d2mdtdf(i)=0.d0
            d2mdt2(i)=0.d0
            d2mdf2(i)=0.d0
            totm(i)=0.d0
            do lm=1,ndimv
              if(symb(lm,i,isort)) then
                isymb=lmsymb(lm,i,isort)
                mtb=indmtb(isymb,isort)
                ro0=spmt(mtb+irad)
                ro1=dspdr(irad,isymb)
                ro2=d2spdr2(irad,isymb)
                totm(i)=totm(i)+ro0*yl(lm)
                dmdr(i)=dmdr(i)+ro1*yl(lm)
                dmdt(i)=dmdt(i)+ro0*d12yl(lm,1)
                dmdf(i)=dmdf(i)+ro0*d12yl(lm,2)
                d2mdr2(i)=d2mdr2(i)+ro2*yl(lm)
                d2mdt2(i)=d2mdt2(i)+ro0*d12yl(lm,3)
                d2mdf2(i)=d2mdf2(i)+ro0*d12yl(lm,5)
                d2mdrdt(i)=d2mdrdt(i)+ro1*d12yl(lm,1)
                d2mdrdf(i)=d2mdrdf(i)+ro1*d12yl(lm,2)
                d2mdtdf(i)=d2mdtdf(i)+ro0*d12yl(lm,4)
              endif
            enddo  !!! over lm
          enddo
          sm=dot_product(totm,b_ext)
          if(sm.lt.1.d-12) then
            dsmdr=0.d0
            dsmdt=0.d0
            dsmdf=0.d0
            d2smdr2=0.d0
            d2smdrdt=0.d0
            d2smdrdf=0.d0
            d2smdt2=0.d0
            d2smdf2=0.d0
            d2smdtdf=0.d0
          else
            dsmdr=dot_product(dmdr,b_ext)
            dsmdt=dot_product(dmdt,b_ext)
            dsmdf=dot_product(dmdf,b_ext)
            d2smdr2=dot_product(d2mdr2,b_ext)
            d2smdrdt=dot_product(d2mdrdt,b_ext)
            d2smdrdf=dot_product(d2mdrdf,b_ext)
            d2smdt2=dot_product(d2mdt2,b_ext)
            d2smdf2=dot_product(d2mdf2,b_ext)
            d2smdtdf=dot_product(d2mdtdf,b_ext)
          endif
          sm_lap=2.d0*r1*dsmdr+d2smdr2+costet*r1*rst*dsmdt+r2*d2smdt2
     &               +rst2*d2smdf2
c ------- Forming RO_up and RO_down values ----------------------------
          rotot(2)=rotot(1)+0.5d0*sm
          rotot(1)=rotot(1)-0.5d0*sm
          rolap(2)=rolap(1)+0.5d0*sm_lap
          rolap(1)=rolap(1)-0.5d0*sm_lap
          ror(2)=ror(1)+0.5d0*dsmdr
          ror(1)=ror(1)-0.5d0*dsmdr
          rot(2)=rot(1)+0.5d0*dsmdt
          rot(1)=rot(1)-0.5d0*dsmdt
          rof(2)=rof(1)+0.5d0*dsmdf
          rof(1)=rof(1)-0.5d0*dsmdf
          rorr(2)=rorr(1)+0.5d0*d2smdr2
          rorr(1)=rorr(1)-0.5d0*d2smdr2
          rort(2)=rort(1)+0.5d0*d2smdrdt
          rort(1)=rort(1)-0.5d0*d2smdrdt
          rorf(2)=rorf(1)+0.5d0*d2smdrdf
          rorf(1)=rorf(1)-0.5d0*d2smdrdf
          rott(2)=rott(1)+0.5d0*d2smdt2
          rott(1)=rott(1)-0.5d0*d2smdt2
          rotf(2)=rotf(1)+0.5d0*d2smdtdf
          rotf(1)=rotf(1)-0.5d0*d2smdtdf
          roff(2)=roff(1)+0.5d0*d2smdf2
          roff(1)=roff(1)-0.5d0*d2smdf2
        endif
      endif
      do ispin=1,2
        gro(1,ispin)=ror(ispin)
        gro(2,ispin)=r1*rot(ispin)
        gro(3,ispin)=rst*rof(ispin)
        g(1)=sintet*cosphi*gro(1,ispin)+costet*cosphi*gro(2,ispin)
     &                                 -sinphi*gro(3,ispin)
        g(2)=sintet*sinphi*gro(1,ispin)+costet*sinphi*gro(2,ispin)
     &                                 +cosphi*gro(3,ispin)
        g(3)=costet*gro(1,ispin)-sintet*gro(2,ispin)
        gro(:,ispin)=g
        romod(ispin)=ddot(3,gro(1,ispin),1,gro(1,ispin),1)
        romod(ispin)=dsqrt(romod(ispin))
        gmod(1,ispin)=(ror(ispin)*rorr(ispin)+r2*rot(ispin)*
     &      (rort(ispin)-r1*rot(ispin))+rst2*rof(ispin)*
     &      (rorf(ispin)-r1*rof(ispin)))/romod(ispin)
        gmod(2,ispin)=(ror(ispin)*rort(ispin)+r2*rott(ispin)*
     &      rot(ispin)+rst2*rof(ispin)*(rotf(ispin)-
     &      costet/sintet*rof(ispin)))/romod(ispin)
        gmod(3,ispin)=(ror(ispin)*rorf(ispin)+r2*rotf(ispin)
     &         *rot(ispin)+rst2*rof(ispin)*roff(ispin))/romod(ispin)
        gmod(2,ispin)=r1*gmod(2,ispin)
        gmod(3,ispin)=rst*gmod(3,ispin)
      enddo
      do i=1,3
        gro(i,3)=gro(i,1)+gro(i,2)
      enddo
      romod(3)=0.d0
      do i=1,3
        romod(3)=romod(3)+gro(i,3)**2
      enddo
      romod(3)=dsqrt(romod(3))
      rrr=rorr(1)+rorr(2)
      rrt=rort(1)+rort(2)
      rrf=rorf(1)+rorf(2)
      rtt=rott(1)+rott(2)
      rtf=rotf(1)+rotf(2)
      rff=roff(1)+roff(2)
      rvv=ror(1)+ror(2)
      rt=rot(1)+rot(2)
      rf=rof(1)+rof(2)
      gmod(1,3)=(rvv*rrr+r2*rt*(rrt-r1*rt)+rst2*rf*(rrf-r1*rf))
     &                                        /romod(3)
      gmod(2,3)=(rvv*rrt+r2*rtt*rt+rst2*rf*(rtf-costet/sintet*rf))
     &                                                    /romod(3)
      gmod(3,3)=(rvv*rrf+r2*rtf*rt+rst2*rf*rff)/romod(3)
      gmod(2,3)=r1*gmod(2,3)
      gmod(3,3)=rst*gmod(3,3)
c ----- Transform Grad |Grad Rho| from Spherical to Cartesian ------
      do ispin=1,3
        g(1)=sintet*cosphi*gmod(1,ispin)+costet*cosphi*gmod(2,ispin)
     &                                  -sinphi*gmod(3,ispin)
        g(2)=sintet*sinphi*gmod(1,ispin)+costet*sinphi*gmod(2,ispin)
     &                                  +cosphi*gmod(3,ispin)
        g(3)=costet*gmod(1,ispin)-sintet*gmod(2,ispin)
        gmod(:,ispin)=g
      enddo
      end
