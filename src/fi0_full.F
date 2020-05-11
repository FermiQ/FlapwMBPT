      subroutine fi0_full(fifim,isp1,isp2)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: isp1,isp2
      real*8, intent(out) :: fifim(maxel,maxel,n_pbmtm,nsort)
      integer :: isort,ndimb,ndimpb,lm3,l3,nm3,lm,l,nm,ii,iii,indx,lget,
     &	         icg1,icg2,icg,lmi,l_i,iei,lfi,ie3,km3,lf3,ie,km,lf,
     &           in3,in,jj3,jj,i,li3,ii3,mj3,li,mj,ispi1,ispi2
      real*8 :: cgg
      ispi1=min(isp1,nspin_1)
      ispi2=min(isp2,nspin_1)
      fifim=0.d0
      do isort=1,nsort
        ndimb=nrel*(lmb(isort)+1)**2
        ndimpb=(lmpb(isort)+1)**2
        do lm3=1,ndimb
          if(irel.le.1) then
            l3=lget(lm3)
            li3=l3+1
          else if(irel.eq.2) then
            call getlimj(lm3,l3,ii3,mj3,li3,0)
          endif
          nm3=ntle(l3,isort)
          do lm=1,ndimb
            if(irel.le.1) then
              l=lget(lm)
              li=l+1
            else if(irel.eq.2) then
              call getlimj(lm,l,ii,mj,li,0)
            endif
            nm=ntle(l,isort)
            if(irel.le.1) then
              ii = max0(lm3,lm)
              iii = min0(lm3,lm)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
            else
              icg1=indxcgr(lm3,lm,1,0)
              icg2=indxcgr(lm3,lm,2,0)
            endif
            do icg = icg1, icg2
              if(irel.le.1) then
                lmi = jcg(icg)
                cgg=cg(icg)
              else
                lmi=jcgr(icg)
                cgg=cgr(icg)
              endif
              if(lmi.le.ndimpb) then
                l_i=lget(lmi)
                do iei=1,ntle_pb(l_i,isort)
                  lfi=ind_prod0(iei,l_i,isort)
                  i=ind_prod(iei,lmi,isort)
                  do ie3=1,nm3
                    in3=1
                    if(augm(ie3,l3,isort)=='LAW') in3=2
                    if(augm(ie3,l3,isort)=='AWL') in3=2
                    do jj3=1,in3
                      km3=indbasa(jj3,ie3,lm3,isort)
                      lf3=lf_isz(km3,isort)
                      do ie=1,nm
                        in=1
                        if(augm(ie,l,isort)=='LAW') in=2
                        if(augm(ie,l,isort)=='AWL') in=2
                        do jj=1,in
                          km=indbasa(jj,ie,lm,isort)
                          lf=lf_isz(km,isort)
                          fifim(km,km3,i,isort)=
     &                          cgg*fi0(lf3,lf,lfi,isort,ispi1,ispi2)
                        enddo
                      enddo   !! over ie
                    enddo
                  enddo   !! over ie3
                enddo   !!! over iei
              endif
            enddo   !! over icg
          enddo   !! over lm
        enddo  !! over lm3
      enddo  !! over isort
      end
      
      subroutine fi0_full_red(fifim,isp1,isp2)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: isp1,isp2
      real*8, intent(out) :: fifim(nrel,maxel_red,maxel_red,n_pbmtm_red,
     &                             nsort)
      integer :: isort,ndimb,ndimpb,lm3,l3,nm3,lm,l,nm,ii,iii,indx,lget,
     &	       icg1,icg2,icg,lmi,l_i,iei,lfi,ie3,km3,lf3,ie,km,lf,
     &           in3,in,jj3,jj,i,li3,ii3,mj3,li,mj,ispi1,ispi2
      real*8 :: cgg,cggi
      ispi1=min(isp1,nspin_1)
      ispi2=min(isp2,nspin_1)
      fifim=0.d0
      do isort=1,nsort
        ndimb=nrel*(lmb_red(isort)+1)**2
        ndimpb=(lmpb_red(isort)+1)**2
        do lm3=1,ndimb
          if(irel.le.1) then
            l3=lget(lm3)
            li3=l3+1
          else if(irel.eq.2) then
            call getlimj(lm3,l3,ii3,mj3,li3,0)
          endif
          nm3=ntle(l3,isort)
          do lm=1,ndimb
            if(irel.le.1) then
              l=lget(lm)
              li=l+1
            else if(irel.eq.2) then
              call getlimj(lm,l,ii,mj,li,0)
            endif
            nm=ntle(l,isort)
            if(irel.le.1) then
              ii = max0(lm3,lm)
              iii = min0(lm3,lm)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
            else
              icg1=indxcgr(lm3,lm,1,0)
              icg2=indxcgr(lm3,lm,2,0)
            endif
            do icg = icg1, icg2
              if(irel.le.1) then
                lmi = jcg(icg)
                cgg=cg(icg)
              else
                lmi=jcgr(icg)
                cgg=real(cgr(icg))
                cggi=imag(cgr(icg))
              endif
              if(lmi.le.ndimpb) then
                l_i=lget(lmi)
                do iei=1,ntle_pb_red(l_i,isort)
                  lfi=ind_prod0_red(iei,l_i,isort)
                  i=ind_prod_red(iei,lmi,isort)
                  do ie3=1,nm3
                    if(correlated(ie3,l3,isort)=='N') cycle
                    in3=1
                    if(augm(ie3,l3,isort)=='LAW') in3=2
                    if(augm(ie3,l3,isort)=='AWL') in3=2
                    if(correlated(ie3,l3,isort)=='A') in3=1
                    do jj3=1,in3
                      km3=indbasa_red(jj3,ie3,lm3,isort)
                      lf3=lf_isz_red(km3,isort)
                      do ie=1,nm
                        if(correlated(ie,l,isort)=='N') cycle
                        in=1
                        if(augm(ie,l,isort)=='LAW') in=2
                        if(augm(ie,l,isort)=='AWL') in=2
                        if(correlated(ie,l,isort)=='A') in=1
                        do jj=1,in
                          km=indbasa_red(jj,ie,lm,isort)
                          lf=lf_isz_red(km,isort)
                          fifim(1,km,km3,i,isort)=
     &                        cgg*fi0_red(lf3,lf,lfi,isort,ispi1,ispi2)
                          if(irel==2) fifim(2,km,km3,i,isort)=
     &                        cggi*fi0_red(lf3,lf,lfi,isort,ispi1,ispi2)
                        enddo
                      enddo   !! over ie
                    enddo
                  enddo   !! over ie3
                enddo   !!! over iei
              endif
            enddo   !! over icg
          enddo   !! over lm
        enddo  !! over lm3
      enddo  !! over isort
      end
      
      subroutine fi0_full_red_1(fifim,isp1,isp2)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: isp1,isp2
      real*8, intent(out) :: fifim(nrel,maxel,maxel,n_pbmtm_red,nsort)
      integer :: isort,ndimb,ndimpb,lm3,l3,nm3,lm,l,nm,ii,iii,indx,lget,
     &	       icg1,icg2,icg,lmi,l_i,iei,lfi,ie3,km3,lf3,ie,km,lf,
     &           in3,in,jj3,jj,i,li3,ii3,mj3,li,mj,ispi1,ispi2
      real*8 :: cgg,cggi
      ispi1=min(isp1,nspin_1)
      ispi2=min(isp2,nspin_1)
      fifim=0.d0
      do isort=1,nsort
        ndimb=nrel*(lmb(isort)+1)**2
        ndimpb=(lmpb_red(isort)+1)**2
        do lm3=1,ndimb
          if(irel.le.1) then
            l3=lget(lm3)
            li3=l3+1
          else if(irel.eq.2) then
            call getlimj(lm3,l3,ii3,mj3,li3,0)
          endif
          nm3=ntle(l3,isort)
          do lm=1,ndimb
            if(irel.le.1) then
              l=lget(lm)
              li=l+1
            else if(irel.eq.2) then
              call getlimj(lm,l,ii,mj,li,0)
            endif
            nm=ntle(l,isort)
            if(irel.le.1) then
              ii = max0(lm3,lm)
              iii = min0(lm3,lm)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
            else
              icg1=indxcgr(lm3,lm,1,0)
              icg2=indxcgr(lm3,lm,2,0)
            endif
            do icg = icg1, icg2
              if(irel.le.1) then
                lmi = jcg(icg)
                cgg=cg(icg)
              else
                lmi=jcgr(icg)
                cgg=real(cgr(icg))
                cggi=imag(cgr(icg))
              endif
              if(lmi.le.ndimpb) then
                l_i=lget(lmi)
                do iei=1,ntle_pb_red(l_i,isort)
                  lfi=ind_prod0_red(iei,l_i,isort)
                  i=ind_prod_red(iei,lmi,isort)
                  do ie3=1,nm3
                    in3=1
                    if(augm(ie3,l3,isort)=='LAW') in3=2
                    if(augm(ie3,l3,isort)=='AWL') in3=2
                    do jj3=1,in3
                      km3=indbasa(jj3,ie3,lm3,isort)
                      lf3=lf_isz(km3,isort)
                      do ie=1,nm
                        in=1
                        if(augm(ie,l,isort)=='LAW') in=2
                        if(augm(ie,l,isort)=='AWL') in=2
                        do jj=1,in
                          km=indbasa(jj,ie,lm,isort)
                          lf=lf_isz(km,isort)
                          fifim(1,km,km3,i,isort)=
     &                      cgg*fi0_f_red(lf3,lf,lfi,isort,ispi1,ispi2)
                          if(irel==2) fifim(2,km,km3,i,isort)=
     &                      cggi*fi0_f_red(lf3,lf,lfi,isort,ispi1,ispi2)
                        enddo
                      enddo   !! over ie
                    enddo
                  enddo   !! over ie3
                enddo   !!! over iei
              endif
            enddo   !! over icg
          enddo   !! over lm
        enddo  !! over lm3
      enddo  !! over isort
      end
