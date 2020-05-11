      subroutine sym_sig_rs_c_rel(f_red,f)
c     Get F for all R-sites from symmetrized (reduced) representation
c     i.e. to get F(R,tEL,t'E'L') from F_RED(EL,E'L',Rtt')
c     where Rtt' is reduced to the compact number of terms
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(in) :: f_red(maxel_red,maxel_red,nrr_max_c)
      complex*16, intent(out) :: f(nfun_red,nfun_red,nqdiv_c)
      integer :: j,ig,iatom,isort,jatom,jsort,jnd0,ind0,ir,ir0,iatom0,
     &           jatom0,l,ie,lm,km,lm0,i,in,jn,li,ist,ii,jj,mj
      complex*16 :: y0j(2*maxb+2),t0j(2*maxb+2)
      complex*16, allocatable :: tmpj(:,:)
      allocate(tmpj(maxel_red,maxel_red))
      f=(0.d0,0.d0)
      do jatom=1,natom
        jsort=is(jatom)
        jnd0=io_lem_red(jatom)-1
        do iatom=1,natom
          isort=is(iatom)
          ind0=io_lem_red(iatom)-1
          do ir=1,nqdiv_c
            ir0=r0_pair_c(ir,iatom,jatom)
            iatom0=nrr_red_c(2,ir0)
            jatom0=nrr_red_c(3,ir0)
            ig=ig_pair_c(ir,iatom,jatom)
            do j=1,lfunm_red(jsort)
              li=0
              do l=0,lmb_red(isort)
                ist=-1
                if(l.eq.0) ist=1
                do ii=ist,1,2
                  li=li+1
                  jj=l+l+ii
                  do ie=1,ntle(l,isort)
                    if(correlated(ie,l,isort)=='N') cycle
                    in=1
                    if(augm(ie,l,isort)=='LAW') in=2
                    if(augm(ie,l,isort)=='AWL') in=2
                    if(correlated(ie,l,isort)=='A') in=1
                    do jn=1,in
                      lm0=0
                      do mj=-jj,jj,2
                        lm0=lm0+1
                        call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa_red(jn,ie,lm,isort)
                        t0j(lm0)=f_red(km,j,ir0)
                      enddo
                      call rotate2(t0j,y0j,l,ii,uj(1,ig),1,1)
                      do mj=-jj,jj,2
                        call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa_red(jn,ie,lm,isort)
                        lm0=(jj+mj)/2+1
                        tmpj(km,j)=y0j(lm0)
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo   !! over j
            do i=1,lfunm_red(isort)
              li=0
              do l=0,lmb_red(jsort)
                ist=-1
                if(l.eq.0) ist=1
                do ii=ist,1,2
                  li=li+1
                  jj=l+l+ii
                  do ie=1,ntle(l,jsort)
                    if(correlated(ie,l,jsort)=='N') cycle
                    in=1
                    if(augm(ie,l,jsort)=='LAW') in=2
                    if(augm(ie,l,jsort)=='AWL') in=2
                    if(correlated(ie,l,jsort)=='A') in=1
                    do jn=1,in
                      lm0=0
                      do mj=-jj,jj,2
                        lm0=lm0+1
                        call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa_red(jn,ie,lm,jsort)
                        t0j(lm0)=tmpj(i,km)
                      enddo
                      call rotate2(t0j,y0j,l,ii,uj(1,ig),1,2)
                      do mj=-jj,jj,2
                        call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa_red(jn,ie,lm,jsort)
                        lm0=(jj+mj)/2+1
                        f(ind0+i,jnd0+km,ir)=
     &                      f(ind0+i,jnd0+km,ir)+y0j(lm0)
                      enddo
                    enddo
                  enddo   !! over ie
                enddo   !! over ii
              enddo   !! over l
            enddo  !! over i
          enddo
        enddo
      enddo
      deallocate(tmpj)
      end
