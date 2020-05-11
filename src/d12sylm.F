      subroutine d12sylm(sintet,costet,sinfi,cosfi,d12yl,lmax,yl)
      use atom_mod
      implicit none
      integer :: lmax,l,ind0,m,lm,ind1
      real*8 :: yl(limlpb),d12yl(limlpb,5),sintet,costet,
     &        sinfi,cosfi,s2tet,ctgtet,sl,s
      s2tet=sintet**2
      ctgtet=costet/sintet
      do l=0,lmax
         ind0=l*(l+1)+1
         if(l.gt.0) sl=-dsqrt((ind0-1)/2.d0)
         do m=-l,l
            if(m.eq.0) then
               if(l.eq.0) d12yl(1,1)=0.d0
               if(l.gt.0) d12yl(ind0,1)=sl*(cosfi*yl(ind0+1)+
     &                                      sinfi*yl(ind0-1))
            else if(m.lt.0) then
               d12yl(ind0+m,1)=-m*ctgtet*yl(ind0+m)
               lm=(l+m)*(l-m+1)
               if(lm.ne.0) d12yl(ind0+m,1)=d12yl(ind0+m,1)+
     &            dsqrt(dfloat(lm))*(sinfi*yl(ind0-m+1)-
     &                               cosfi*yl(ind0+m-1))
            else
               d12yl(ind0+m,1)=m*ctgtet*yl(ind0+m)
               lm=(l-m)*(l+m+1)
               if(lm.ne.0) d12yl(ind0+m,1)=d12yl(ind0+m,1)-
     &            dsqrt(dfloat(lm))*(sinfi*yl(ind0-m-1)+
     &                               cosfi*yl(ind0+m+1))
            endif
            d12yl(ind0+m,2)=-m*yl(ind0-m)
            d12yl(ind0+m,5)=-m*m*yl(ind0+m)
         enddo  !!! over m
      enddo     !!! over l
      do l=0,lmax
         ind0=l*(l+1)+1
         if(l.gt.0) sl=-dsqrt((ind0-1)/2.d0)
         do m=-l,-1
            ind1=ind0+m
            d12yl(ind1,3)=m*(yl(ind1)/s2tet-ctgtet*d12yl(ind1,1))
            d12yl(ind1,4)=m**2*ctgtet*yl(ind0-m)
            lm=(l+m)*(l-m+1)
            if(lm.ne.0) then
               s=dsqrt(dfloat(lm))
               d12yl(ind1,3)=d12yl(ind1,3)+s*(sinfi*d12yl(ind0-m+1,1)-
     &            cosfi*d12yl(ind0+m-1,1))
               d12yl(ind1,4)=d12yl(ind1,4)+s*m*(sinfi*yl(ind0+m-1)+
     &            cosfi*yl(ind0-m+1))
            endif
         enddo
         if(l.eq.0) d12yl(ind0,3)=0.d0
         if(l.gt.0) d12yl(ind0,3)=sl*(cosfi*d12yl(ind0+1,1)+sinfi*
     &                 d12yl(ind0-1,1))
         d12yl(ind0,4)=0.d0
         do m=1,l
            ind1=ind0+m
            d12yl(ind1,3)=-m*(yl(ind1)/s2tet-ctgtet*d12yl(ind1,1))
            d12yl(ind1,4)=-m**2*ctgtet*yl(ind0-m)
            lm=(l-m)*(l+m+1)
            if(lm.ne.0) then
               s=dsqrt(dfloat(lm))
               d12yl(ind1,3)=d12yl(ind1,3)-s*(sinfi*d12yl(ind0-m-1,1)+
     &            cosfi*d12yl(ind0+m+1,1))
               d12yl(ind1,4)=d12yl(ind1,4)-s*m*(sinfi*yl(ind0+m+1)-
     &            cosfi*yl(ind0-m-1))
            endif
         enddo
      enddo
      end
