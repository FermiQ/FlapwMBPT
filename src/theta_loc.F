      subroutine theta_loc(lda,roup,rodn,gro,gmod,romod,rolap,th,key)
c     key = 0 - result is provided in bosonic 4-vector
c     key = 1 - the result is provided in spin-spin 2*2 form
	  use atom_mod
	  use manager_mod
	  use solid_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: lda,key
      real*8, intent(in) :: roup,rodn,gro(3,3),gmod(3,3),romod(3),
     &                      rolap(2)
      real*8, intent(out) :: th(n_ixc_0)
      integer :: i,ii
	  real*8 :: a1,a2,b1,b2,excup,excdn,fac,hh,hh2,thet(2,2),a1t,a2t,
     &          b1t,b2t
      thet=0.d0
      fac=1.d-3
      if(roup>1.d-8) then
        hh=roup*fac
        hh2=hh+hh
        call exch_corr(roup+hh,rodn,lda,nloc,gro,gmod,romod,rolap,EXCUP,
     &                 EXCDN,b1,b2)
        call exch_corr(roup-hh,rodn,lda,nloc,gro,gmod,romod,rolap,EXCUP,
     &                 EXCDN,a1,a2)
        if(key_hybrid==1) then
          if(hybrid_type<=2) then
            call exch_corr(roup+hh,rodn,6,nloc,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,b1t,b2t)
          else
            call exch_corr(roup+hh,rodn,6,3,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,b1t,b2t)
          endif
          b1=b1-alpha_hbr*b1t
          b2=b2-alpha_hbr*b2t
          if(hybrid_type<=2) then
            call exch_corr(roup-hh,rodn,6,nloc,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,a1t,a2t)
          else
            call exch_corr(roup-hh,rodn,6,3,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,a1t,a2t)
          endif
          a1=a1-alpha_hbr*a1t
          a2=a2-alpha_hbr*a2t
        endif
        thet(2,2)=(b1-a1)/hh2
        thet(1,2)=(b2-a2)/hh2
      endif
      if(rodn>1.d-8) then
        hh=rodn*fac
        hh2=hh+hh
        call exch_corr(roup,rodn+hh,lda,nloc,gro,gmod,romod,rolap,EXCUP,
     &                 EXCDN,b1,b2)
        call exch_corr(roup,rodn-hh,lda,nloc,gro,gmod,romod,rolap,EXCUP,
     &                 EXCDN,a1,a2)
        if(key_hybrid==1) then
          if(hybrid_type<=2) then
            call exch_corr(roup,rodn+hh,6,nloc,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,b1t,b2t)
          else
            call exch_corr(roup,rodn+hh,6,3,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,b1t,b2t)
          endif
          b1=b1-alpha_hbr*b1t
          b2=b2-alpha_hbr*b2t
          if(hybrid_type<=2) then
            call exch_corr(roup,rodn-hh,6,nloc,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,a1t,a2t)
          else
            call exch_corr(roup,rodn-hh,6,3,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,a1t,a2t)
          endif
          a1=a1-alpha_hbr*a1t
          a2=a2-alpha_hbr*a2t
        endif
        thet(2,1)=(b1-a1)/hh2
        thet(1,1)=(b2-a2)/hh2
      endif
      if(key==1) then
        th(1)=thet(1,1)
        th(2)=thet(2,1)
        th(3)=thet(1,2)
        th(4)=thet(2,2)
      else if(key==0) then
        do i=1,n_ixc_0
          ii=ind_ixc(i)
          if(ii==1) then
            th(i)=0.25d0*(thet(1,1)+thet(2,2)+thet(1,2)+thet(2,1))
          else if(ii==2) then
            th(i)=0.25d0*(thet(1,1)+thet(2,2)-thet(1,2)-thet(2,1))
          else if(ii==3) then
            th(i)=0.25d0*(thet(2,2)-thet(1,1))
          else if(ii==4) then
            call exch_corr(roup,rodn,lda,nloc,gro,gmod,romod,rolap,
     &                     EXCUP,EXCDN,b1,b2)
            if(key_hybrid==1) then
              if(hybrid_type<=2) then
                call exch_corr(roup,rodn,6,nloc,gro,gmod,romod,rolap,
     &                         EXCUP,EXCDN,b1t,b2t)
              else
                call exch_corr(roup,rodn,6,3,gro,gmod,romod,rolap,
     &                         EXCUP,EXCDN,b1t,b2t)
              endif
              b1=b1-alpha_hbr*b1t
              b2=b2-alpha_hbr*b2t
            endif
            if(abs(roup-rodn)<1.d-8) then
              th(i)=0.d0
            else
              th(i)=0.5d0*(b1-b2)/(roup-rodn) ! Theta{+-}=1/2 Theta^{XX}
            endif
          endif
        enddo
      endif
      end
