      subroutine c_alpha_int(nb,x,k)
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nb,k
      complex*16, intent(out) :: x(nb,nb,3)
      integer :: jgb,jgb0,igb,igb0,ia,ib,ic,ind,jp,jh,ipp,ih,nbndt
      real*8 :: gki(3),gkj(3),pi2a,c1,c05,c2,g2j,enj,c2j,anormj,
     &          g2i,eni,c2i,anormi,ann,pi2
      complex*16 :: ss,const1,const2
      pi2=pi+pi
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      x=(0.d0,0.d0)
      pi2a=pi2/par
      c1=amegaint/amega
      nbndt=nb/2
      c05=0.5d0*clight
      c2=clight*clight
      do jgb=1,nbndt
        jgb0=indgb(jgb,k)
        gkj=pi2a*(pnt(:,k)+gbs(:,jgb0))
        g2j=dot_product(gkj,gkj)
        enj=c05*(sqrt(c2+4.d0*g2j)-clight)
        c2j=c2+enj
        anormj=c2j/dsqrt(c2j*c2j+c2*g2j)
        do igb=1,nbndt
          igb0=indgb(igb,k)
          gki=pi2a*(pnt(:,k)+gbs(:,igb0))
          g2i=dot_product(gki,gki)
          eni=c05*(sqrt(c2+4.d0*g2i)-clight)
          c2i=c2+eni
          anormi=c2i/sqrt(c2i*c2i+c2*g2i)
          ia=igbs(1,jgb0)-igbs(1,igb0)
          ib=igbs(2,jgb0)-igbs(2,igb0)
          ic=igbs(3,jgb0)-igbs(3,igb0)
          ind=indplw(ia,ib,ic)
          if(complex_ro) then
            ss=dcmplx(sovr(ind),sovi(ind))/amega
          else
            ss=sovr(ind)/amega
          endif
          ann=anormi*anormj
          const1=ann*ss*c2/c2i
          const2=ann*ss*c2/c2j
          do jp=1,2
            jh=(jp-1)*nbndt+jgb
            do ipp=1,2
              ih=(ipp-1)*nbndt+igb
c     &&&&& Begin  IV-40-24 &&&&&&&&&&&&&&&&&&&
              if(ipp.eq.1.and.jp.eq.1) then
                x(ih,jh,1)=const2*dcmplx(gkj(1),-gkj(2))
     &                    +const1*dcmplx(gki(1),gki(2))
                x(ih,jh,2)=const2*dcmplx(gkj(2),gkj(1))
     &                    +const1*dcmplx(gki(2),-gki(1))
                x(ih,jh,3)=const2*gkj(3)+const1*gki(3)
              else if(ipp.eq.1.and.jp.eq.2) then
                x(ih,jh,1)=const2*dcmplx(gkj(3),0.d0)
     &                    +const1*dcmplx(-gki(3),0.d0)
                x(ih,jh,2)=const2*dcmplx(0.d0,gkj(3))
     &                    +const1*dcmplx(0.d0,-gki(3))
                x(ih,jh,3)=const2*dcmplx(-gkj(1),-gkj(2))
     &                    +const1*dcmplx(gki(1),gkj(2))
              else if(ipp.eq.2.and.jp.eq.1) then
                x(ih,jh,1)=const2*dcmplx(-gkj(3),0.d0)
     &                    +const1*dcmplx(gki(3),0.d0)
                x(ih,jh,2)=const2*dcmplx(0.d0,gkj(3))
     &                    +const1*dcmplx(0.d0,-gki(3))
                x(ih,jh,3)=const2*dcmplx(gkj(1),-gkj(2))
     &                    +const1*dcmplx(-gki(1),gkj(2))
              else if(ipp.eq.2.and.jp.eq.2) then
                x(ih,jh,1)=const2*dcmplx(gkj(1),gkj(2))
     &                    +const1*dcmplx(gki(1),-gki(2))
                x(ih,jh,2)=const2*dcmplx(gkj(2),-gkj(1))
     &                    +const1*dcmplx(gki(2),gki(1))
                x(ih,jh,3)=const2*gkj(3)+const1*gki(3)
              endif
            enddo  !!! over ip
          enddo  !!! over jp
        enddo  !!! over igb
      enddo  !!! over jgb
      end
