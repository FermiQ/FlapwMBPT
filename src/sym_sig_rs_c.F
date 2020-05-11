      subroutine sym_sig_rs_c(f_red,f)
c     Get F for all R-sites from symmetrized (reduced) representation
c     i.e. to get F(R,tEL,t'E'L') from F_RED(EL,E'L',Rtt')
c     where Rtt' is reduced to compact number of terms
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(in) :: f_red(maxel_red,maxel_red,nrr_max_c)
      complex*16, intent(out) :: f(nfun_red,nfun_red,nqdiv_c)
      integer :: j,ig,iatom,isort,jatom,jsort,jnd0,ind0,ir,ir0,l,lp1,ie,
     &           m,lm,km,lm0,i,in,jn,it
      real*8 :: y0(2*maxb+1),t0(2*maxb+1)
      real*8, allocatable :: tmp(:,:,:)
      allocate(tmp(maxel_red,maxel_red,2))
      f=(0.d0,0.d0)
      do jatom=1,natom
        jsort=is(jatom)
        jnd0=io_lem_red(jatom)-1
        do iatom=1,natom
          isort=is(iatom)
          ind0=io_lem_red(iatom)-1
          do ir=1,nqdiv_c
            ir0=r0_pair_c(ir,iatom,jatom)
            ig=ig_pair_c(ir,iatom,jatom)
            do j=1,lfunm_red(jsort)
              do l=0,lmb_red(isort)
                lp1=l+1
                do ie=1,ntle(l,isort)
                  if(correlated(ie,l,isort)=='N') cycle
                  in=1
                  if(augm(ie,l,isort)=='LAW') in=2
                  if(augm(ie,l,isort)=='AWL') in=2
                  if(correlated(ie,l,isort)=='A') in=1
                  do jn=1,in
                    do it=1,2  ! real and imag parts
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=indbasa_red(jn,ie,lm,isort)
                        lm0=l+m+1
                        if(it==1) t0(lm0)=real(f_red(km,j,ir0))
                        if(it==2) t0(lm0)=imag(f_red(km,j,ir0))
                      enddo
                      call rotate1(t0,y0,l,u(1,ig),1)
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=indbasa_red(jn,ie,lm,isort)
                        lm0=l+m+1
                        tmp(km,j,it)=y0(lm0)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over l
            enddo   !! over j
            do i=1,lfunm_red(isort)
              do l=0,lmb_red(jsort)
                lp1=l+1
                do ie=1,ntle(l,jsort)
                  if(correlated(ie,l,jsort)=='N') cycle
                  in=1
                  if(augm(ie,l,jsort)=='LAW') in=2
                  if(augm(ie,l,jsort)=='AWL') in=2
                  if(correlated(ie,l,jsort)=='A') in=1
                  do jn=1,in
                    do it=1,2  ! real and imag parts
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=indbasa_red(jn,ie,lm,jsort)
                        lm0=l+m+1
                        t0(lm0)=tmp(i,km,it)
                      enddo
                      call rotate1(t0,y0,l,u(1,ig),1)
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=indbasa_red(jn,ie,lm,jsort)
                        lm0=l+m+1
                        if(it==1) f(ind0+i,jnd0+km,ir)=
     &                    f(ind0+i,jnd0+km,ir)+y0(lm0)
                        if(it==2) f(ind0+i,jnd0+km,ir)=
     &                    f(ind0+i,jnd0+km,ir)+(0.d0,1.d0)*y0(lm0)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over l
            enddo  !! over i
          enddo
        enddo
      enddo
      deallocate(tmp)
      end
