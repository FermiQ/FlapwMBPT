      subroutine radli(v0,isort,ispin,acc,dff,nds)
c     Radial equation solver for different IDMD:
c     IDMD = 0 - for the given energy (center of occupied bands)
c     IDMD = 1 - for the given NODES and Logarithmic derivative
c     IDMD = 2 - Phi(Smt)=0; Phi'(Smt)/=0 - bound sol-n, but not deriv.
c     IDMD = 3 - Phi(Smt)=Phi'(Smt)=0 - totally bound sol-n (for deep E)
c     IDMD = 4 - for the given energy (equal to the Chemical Potential)
c     IDMD = 5 - for the given energy (center of valence bands)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: isort,ispin
      real*8, intent(in) :: v0(0:maxnrad)
      integer, intent(out) :: nds(maxntle,nrel*maxb+1,nspin)
      real*8, intent(out) :: acc(maxntle,nrel*maxb+1,nspin),
     &                       dff(maxntle,nrel*maxb+1,nspin)
      integer :: nres,li,l,i,ntl,ie,nodes,nw,ir,je,nod0,mtw,nod,idm,in,
     &	         jn,in1,jn1,kc,nitec,nre
      real*8 :: sqpi4,c2,e,adet,avv,adv,avs,ads,anorm,fifi,fidfid,dum,
     &	        fid2fid2,dqdall,alf,fll,val,slo,convc,acc2,akap
      real*8, allocatable :: work(:),p0(:),q0(:),p1(:),q1(:),p2(:),
     &                       q2(:),dldot(:),dldot2(:),tmp(:,:,:,:),
     &                       tmp1(:,:,:,:),dtmp(:,:,:,:),
     &                       dp0(:),dq0(:),dp1(:),dq1(:),dp2(:),dq2(:)
      sqpi4=sqrt(4.d0*pi)
      c2=clight*clight
      nw=nrad(isort)
      val=0.d0
      slo=-1.d-10
      allocate(work(0:nw),p0(0:nw),q0(0:nw),p1(0:nw),q1(0:nw),p2(0:nw),
     &         q2(0:nw),tmp(0:nw,2,maxntle,2),tmp1(0:nw,2,maxntle,2),
     &         dldot(maxntle),dldot2(maxntle))
c ------------------------------------------------------------------
      allocate(dp0(0:nw))
      allocate(dq0(0:nw))
      allocate(dp1(0:nw))
      allocate(dq1(0:nw))
      allocate(dp2(0:nw))
      allocate(dq2(0:nw))
      allocate(dtmp(0:nw,2,maxntle,2))
      dff(:,:,ispin)=0.d0
      do li=1,nrel*lmb(isort)+1
        if(irel.le.1) then
          l=li-1
        else
          l=li/2
          i=li-2*l
          if(i.eq.0) i=-1
        endif
        if(irel==1) akap=l*(l+1)
        if(irel==2) then
          if(i.lt.0) akap=l
          if(i.gt.0) akap=-l-1
        endif
        ntl=ntle(l,isort)
        do ie=1,ntl
          nodes=int(ptnl(ie,l,isort,ispin))-l-1           
          idm=idmd(ie,l,isort)
          p0=p_f(0:nw,ie,li,isort,ispin)
          q0=q_f(0:nw,ie,li,isort,ispin)
c -------------- First we solve for Fi -----------------------------         
          e=eny(ie,li,isort,ispin)
          if(iter>1.and.(idm==0.or.idm>=4)) then
c --------- Solve for the given E ------------------------------
            call rad_eq0(e,z(isort),l,i,nod,v0,fifi,r(0,isort),
     &		             nrad(isort),p0,q0,fi(ie,li,isort,ispin),
     &                   dfi(ie,li,isort,ispin),
     &                   dny(ie,li,isort,ispin),dr(0,isort),0,
     &                   h(isort),p1,q1,p1,q1,work)
          else if(idm==3) then
c --------- Solve for the Bound solution (P=0), and (P'=0) ------
            nre=nrad(isort)
            call rad_eq(e+e/2,e-e/2,e,1.d-14,z(isort),l,i,nodes,0.d0,
     &		            -1.d-10,v0,fifi,r(0,isort),nrad(isort),p0,q0,
     &                  dum,fi(ie,li,isort,ispin),
     &                  dfi(ie,li,isort,ispin),kc,nitec,convc,nre,
     &                  dum,dr(0,isort),0,h(isort),
     &                  dny(ie,li,isort,ispin),work,acc2,nod,irel)
            eny(ie,li,isort,ispin)=e
          else if(iter==1.or.idm==1) then
c --------- Solve for the given T ------------------------------ 
            call rad_eq1(e,1.d-14,l,i,ptnl(ie,l,isort,ispin),v0,fifi,
     &		             p0,q0,work,isort,nrad(isort),nod,
     &                   fi(ie,li,isort,ispin),dfi(ie,li,isort,ispin),
     &                   dny(ie,li,isort,ispin))
            eny(ie,li,isort,ispin)=e
          else if(idm==2) then
c --------- Solve for the Bound solution (P=0), but (P'/=0) ------
            call rad_eq4(e,1.d-14,l,i,nodes,v0,fifi,p0,q0,work,isort,
     &                   nrad(isort),nod,fi(ie,li,isort,ispin),
     &                   dfi(ie,li,isort,ispin),dny(ie,li,isort,ispin))
            eny(ie,li,isort,ispin)=e
          endif
          dff(ie,li,ispin)=fifi
          nds(ie,li,ispin)=nod
c --------------------------------------------------------------------------
          q_f(0:nw,ie,li,isort,ispin)=q0
          p_f(0:nw,ie,li,isort,ispin)=p0
          acc(ie,li,ispin)=0.d0
          call dp_dq(e,z(isort),l,r(0,isort),p0,q0,v0,akap,p0,q0,0,nw,
     &               p0,dp0,dq0)
          if(augm(ie,l,isort)/='APW'.and.idm/=3) then
c -------------- Now we solve for Fi_dot ----------------------------- 
            p1=pd_f(0:nw,ie,li,isort,ispin)
            q1=qd_f(0:nw,ie,li,isort,ispin)
            call rad_eq0(e,z(isort),l,i,nod0,v0,fidfid,r(0,isort),
     &	                 nrad(isort),p1,q1,fidot(ie,li,isort,ispin),
     &                   dfidot(ie,li,isort,ispin),dldot(ie),
     &                   dr(0,isort),1,h(isort),p0,q0,p1,q1,work)
            acc(ie,li,ispin)=p1(nw)*q0(nw)-p0(nw)*q1(nw)
            if(irel==1) then
              work(0)=0.d0
              do ir=1,nw
                alf=1.d0+(e+2*z(isort)/r(ir,isort)-v0(ir))/c2
                fll=l*(l+1)/r(ir,isort)**2/alf
                work(ir)=p0(ir)**2*fll/alf*dr(ir,isort)
              enddo
              acc(ie,li,ispin)=acc(ie,li,ispin)
     &		  	            -dqdall(h(isort),work,nrad(isort))/c2
            endif
            pd_f(0:nw,ie,li,isort,ispin)=p1
            qd_f(0:nw,ie,li,isort,ispin)=q1
            call dp_dq(e,z(isort),l,r(0,isort),p1,q1,v0,akap,p0,q0,1,nw,
     &                 p0,dp1,dq1)
c -------------- Now we solve for Fi_dot2 -----------------------------
            if(augm(ie,l,isort)=='LOC') then
              p2=pd2_f(0:nw,ie,li,isort,ispin)
              q2=qd2_f(0:nw,ie,li,isort,ispin)
              call rad_eq0(e,z(isort),l,i,nod0,v0,fid2fid2,r(0,isort),
     &		               nrad(isort),p2,q2,fidot2(ie,li,isort,ispin),
     &                     dfidot2(ie,li,isort,ispin),dldot2(ie),
     &                     dr(0,isort),2,h(isort),p0,q0,p1,q1,work)
              pd2_f(0:nw,ie,li,isort,ispin)=p2
              qd2_f(0:nw,ie,li,isort,ispin)=q2
              call dp_dq(e,z(isort),l,r(0,isort),p2,q2,v0,akap,p1,q1,2,
     &                   nw,p0,dp2,dq2)
            endif
            adet=fi(ie,li,isort,ispin)*dfidot(ie,li,isort,ispin)
     & 	      -dfi(ie,li,isort,ispin)*fidot(ie,li,isort,ispin)
          endif
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          if(augm(ie,l,isort)=='LAW') then
            avv=dfidot(ie,li,isort,ispin)/adet
            adv=-dfi(ie,li,isort,ispin)/adet
            avs=-fidot(ie,li,isort,ispin)/adet
            ads=fi(ie,li,isort,ispin)/adet
            augm_coef(1,li,isort,ispin)=avv
            augm_coef(2,li,isort,ispin)=avs
            augm_coef(3,li,isort,ispin)=0.d0
            do ir=0,nw
              tmp(ir,1,ie,1)=avv*p0(ir)+adv*p1(ir)
              tmp(ir,2,ie,1)=avs*p0(ir)+ads*p1(ir)
              tmp1(ir,1,ie,1)=e*tmp(ir,1,ie,1)+adv*p0(ir)
              tmp1(ir,2,ie,1)=e*tmp(ir,2,ie,1)+ads*p0(ir)
              dtmp(ir,1,ie,1)=avv*dp0(ir)+adv*dp1(ir)
              dtmp(ir,2,ie,1)=avs*dp0(ir)+ads*dp1(ir)
            enddo
            if(irel.ge.1) then
              do ir=0,nrad(isort)
                tmp(ir,1,ie,2)=avv*q0(ir)+adv*q1(ir)
                tmp(ir,2,ie,2)=avs*q0(ir)+ads*q1(ir)
                tmp1(ir,1,ie,2)=e*tmp(ir,1,ie,2)+adv*q0(ir)
                tmp1(ir,2,ie,2)=e*tmp(ir,2,ie,2)+ads*q0(ir)
                dtmp(ir,1,ie,2)=avv*dq0(ir)+adv*dq1(ir)
                dtmp(ir,2,ie,2)=avs*dq0(ir)+ads*dq1(ir)
              enddo
            endif
          else if(augm(ie,l,isort)=='APW') then
            avv=1.d0/fi(ie,li,isort,ispin)
            augm_coef(1,li,isort,ispin)=avv
            augm_coef(2,li,isort,ispin)=0.d0
            augm_coef(3,li,isort,ispin)=0.d0
            do ir=0,nw
              tmp(ir,1,ie,1)=avv*p0(ir)
              tmp1(ir,1,ie,1)=e*tmp(ir,1,ie,1)
              dtmp(ir,1,ie,1)=avv*dp0(ir)
            enddo
            if(irel.ge.1) then
              do ir=0,nrad(isort)
                tmp(ir,1,ie,2)=avv*q0(ir)
                tmp1(ir,1,ie,2)=e*tmp(ir,1,ie,2)
                dtmp(ir,1,ie,2)=avv*dq0(ir)
              enddo
            endif
          else if(augm(ie,l,isort)=='AWL') then
            avv=1.d0/fi(ie,li,isort,ispin)
            avs=fidot(ie,li,isort,ispin)**2
     &         /(fidfid*fi(ie,li,isort,ispin)**2
     &          +fidot(ie,li,isort,ispin)**2)
            avs=sqrt(avs)
            ads=-avs*fi(ie,li,isort,ispin)/fidot(ie,li,isort,ispin)
            augm_coef(1,li,isort,ispin)=avv
            augm_coef(2,li,isort,ispin)=avs
            augm_coef(3,li,isort,ispin)=0.d0
            do ir=0,nw
              tmp(ir,1,ie,1)=avv*p0(ir)
              tmp(ir,2,ie,1)=avs*p0(ir)+ads*p1(ir)
              tmp1(ir,1,ie,1)=e*tmp(ir,1,ie,1)
              tmp1(ir,2,ie,1)=e*tmp(ir,2,ie,1)+ads*p0(ir)
              dtmp(ir,1,ie,1)=avv*dp0(ir)
              dtmp(ir,2,ie,1)=avs*dp0(ir)+ads*dp1(ir)
            enddo
            if(irel.ge.1) then
              do ir=0,nrad(isort)
                tmp(ir,1,ie,2)=avv*q0(ir)
                tmp(ir,2,ie,2)=avs*q0(ir)+ads*q1(ir)
                tmp1(ir,1,ie,2)=e*tmp(ir,1,ie,2)
                tmp1(ir,2,ie,2)=e*tmp(ir,2,ie,2)+ads*q0(ir)
                dtmp(ir,1,ie,2)=avv*dq0(ir)
                dtmp(ir,2,ie,2)=avs*dq0(ir)+ads*dq1(ir)
              enddo
            endif
          else if(augm(ie,l,isort)=='LOC') then
            if(idm==3) then
              do ir=0,nrad(isort)
                tmp(ir,1,ie,1)=p0(ir)
                tmp1(ir,1,ie,1)=e*tmp(ir,1,ie,1)
                dtmp(ir,1,ie,1)=dp0(ir)
              enddo
              if(irel.ge.1) then
                do ir=0,nrad(isort)
                  tmp(ir,1,ie,2)=q0(ir)
                  tmp1(ir,1,ie,2)=e*tmp(ir,1,ie,2)
                  dtmp(ir,1,ie,2)=dq0(ir)
                enddo
              endif
            else
              call local_orb(work,p0,p1,p2,q0,q1,q2,dr(0,isort),
     &	                     nrad(isort),fi(ie,li,isort,ispin),
     &                       dfi(ie,li,isort,ispin),
     &                       fidot(ie,li,isort,ispin),
     &	                     dfidot(ie,li,isort,ispin),
     &                       fidot2(ie,li,isort,ispin),
     &                       dfidot2(ie,li,isort,ispin),adet,
     &                       tmp(0,1,ie,1),tmp(0,1,ie,2),c2,h(isort),
     &                       avv,adv,anorm,fidfid,fid2fid2,
     &                       dtmp(0,1,ie,1),dtmp(0,1,ie,2),
     &                       dp0,dp1,dp2,dq0,dq1,dq2)
              do ir=0,nrad(isort)
                tmp1(ir,1,ie,1)=e*tmp(ir,1,ie,1)
     &		                 +anorm*(adv*p0(ir)+2.d0*p1(ir))
              enddo
              if(irel.ge.1) then
                do ir=0,nrad(isort)
                  tmp1(ir,1,ie,2)=e*tmp(ir,1,ie,2)
     &			               +anorm*(adv*q0(ir)+2.d0*q1(ir))
                enddo
              endif
            endif
          endif
        enddo   !! over ie
c -----------------------------------------------------------------
        do ie=1,ntl
          in=1
          if(augm(ie,l,isort)=='LAW') in=2
          if(augm(ie,l,isort)=='AWL') in=2
          do jn=1,in
            nres=indfun0(jn,ie,li,isort)
            mtw=ind_wf(nres,isort)
            do ir=1,nw
              gfun(mtw+ir,ispin)=tmp(ir,jn,ie,1)/r(ir,isort)
              gfun1(mtw+ir,ispin)=dtmp(ir,jn,ie,1)/r(ir,isort)
     &                           -tmp(ir,jn,ie,1)/r(ir,isort)**2
            enddo
            gfun(mtw,ispin)=gfun(mtw+1,ispin)
            gfun1(mtw,ispin)=gfun1(mtw+1,ispin)
            if(irel.ge.1) then
              do ir=1,nw
                gfund(mtw+ir,ispin)=tmp(ir,jn,ie,2)/r(ir,isort)
                gfund1(mtw+ir,ispin)=dtmp(ir,jn,ie,2)/r(ir,isort)
     &                              -tmp(ir,jn,ie,2)/r(ir,isort)**2
              enddo
              gfund(mtw,ispin)=gfund(mtw+1,ispin)
              gfund1(mtw,ispin)=gfund1(mtw+1,ispin)
            endif
          enddo
        enddo
c ------ Get FFSMT  and  FFHMT ----------------------------------
        do je=1,ntl
          in1=1
          if(augm(je,l,isort)=='LAW') in1=2
          if(augm(je,l,isort)=='AWL') in1=2
          do jn1=1,in1
            do ie=1,ntl
              in=1
              if(augm(ie,l,isort)=='LAW') in=2
              if(augm(ie,l,isort)=='AWL') in=2
              do jn=1,in
                do ir=0,nrad(isort)
                  p0(ir)=tmp(ir,jn,ie,1)*tmp(ir,jn1,je,1)
                  if(irel.ge.1) p0(ir)=p0(ir)+tmp(ir,jn,ie,2)
     &	                                   *tmp(ir,jn1,je,2)/c2
                  p0(ir)=p0(ir)*dr(ir,isort)
                enddo
                ffsmt(jn,jn1,ie,je,li,isort,ispin)=
     &                              dqdall(h(isort),p0,nrad(isort))
                do ir=0,nrad(isort)
                  p0(ir)=tmp(ir,jn,ie,1)*tmp1(ir,jn1,je,1)
                  if(irel.ge.1) p0(ir)=p0(ir)
     &	                          +tmp(ir,jn,ie,2)*tmp1(ir,jn1,je,2)/c2
                  p0(ir)=p0(ir)*dr(ir,isort)
                enddo
                ffhmt(jn,jn1,ie,je,li,isort,ispin)=
     &                dqdall(h(isort),p0,nrad(isort))
              enddo
            enddo
          enddo
        enddo
      enddo   !! over li
      deallocate(dtmp,work,p0,q0,p1,q1,p2,q2,tmp,tmp1,dldot,dldot2,dp0,
     &           dq0,dp1,dq1,dp2,dq2)
      end
