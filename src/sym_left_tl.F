      subroutine sym_left_tl(ig,f0,f,key,li,isort,ns)
c     key=1 conjugated rotation  (D*)     Important only
c     key=2 usual                (D )       for irel=2
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ig,key,isort,ns,li
      complex*16, intent(in) :: f0(ns)
      complex*16, intent(out) :: f(ns)
      integer :: iatom,jatom,ind0,jnd0,l,ie,m,lm,km,lm0,in,jn,iat,ii,jj,
     &           mj
      real*8 :: y0(2*maxb+1),y1(2*maxb+1),t0(2*maxb+1)
      complex*16 :: y0j(2*maxb+2),t0j(2*maxb+2)
      do iat=1,nhsort(isort)
        iatom=iat_sort(iat,isort)
        jatom=ip(iatom,ig)
        ind0=io_lem(iatom)-1   !! zero index for atom iatom
        jnd0=io_lem(jatom)-1
        if(irel<=1) then
          l=li-1
          do ie=1,ntle(l,isort)
            in=1
            if(augm(ie,l,isort)=='LAW') in=2
            if(augm(ie,l,isort)=='AWL') in=2
            do jn=1,in
c ------------- Real components ----------------------------------
              do m=-l,l
                lm=l*(l+1)+m+1
                km=jnd0+indbasa(jn,ie,lm,isort)
                km=nfun_order_back(km)-nfun_ls_adr(li,isort)+1
                lm0=l+m+1
                t0(lm0)=dreal(f0(km))
              enddo
              call rotate1(t0,y0,l,u(1,ig),1)
c ------------- Imaginary components ----------------------------------
              do m=-l,l
                lm=l*(l+1)+m+1
                km=jnd0+indbasa(jn,ie,lm,isort)
                km=nfun_order_back(km)-nfun_ls_adr(li,isort)+1
                lm0=l+m+1
                t0(lm0)=dimag(f0(km))
              enddo
              call rotate1(t0,y1,l,u(1,ig),1)
              do m=-l,l
                lm=l*(l+1)+m+1
                km=ind0+indbasa(jn,ie,lm,isort)
                km=nfun_order_back(km)-nfun_ls_adr(li,isort)+1
                lm0=l+m+1
                f(km)=dcmplx(y0(lm0),y1(lm0))
              enddo
            enddo
          enddo  !! over ie
        else if(irel==2) then
          l=li/2
          ii=2*(li-l-l+2)-5
          jj=l+l+ii
          do ie=1,ntle(l,isort)
            in=1
            if(augm(ie,l,isort)=='LAW') in=2
            if(augm(ie,l,isort)=='AWL') in=2
            do jn=1,in
              lm0=0
              do mj=-jj,jj,2
                lm0=lm0+1
                call getlimj(lm,l,ii,mj,li,1)
                km=jnd0+indbasa(jn,ie,lm,isort)
                km=nfun_order_back(km)-nfun_ls_adr(li,isort)+1
                t0j(lm0)=f0(km)
              enddo
              call rotate2(t0j,y0j,l,ii,uj(1,ig),1,key)
              do mj=-jj,jj,2
                call getlimj(lm,l,ii,mj,li,1)
                km=ind0+indbasa(jn,ie,lm,isort)
                km=nfun_order_back(km)-nfun_ls_adr(li,isort)+1
                lm0=(jj+mj)/2+1
                f(km)=y0j(lm0)
              enddo
            enddo
          enddo
        endif
      enddo   !! over iatom of a given sort
      end
