      subroutine sym_vmi_left_sort(ig,f0,f,n1pbmt,nat1,isort,n1pbm)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ig,n1pbmt,nat1,isort,n1pbm
      complex*16, intent(in) :: f0(n1pbmt)
      complex*16, intent(out) :: f(n1pbmt)
      integer :: iatom,jatom,ind0,jnd0,l,ie,m,lm,km,lm0,iat,jat
      real*8 :: y0(2*maxpb+1),y1(2*maxpb+1),t0(2*maxpb+1)
      do iat=1,nat1
        iatom=iat_sort(iat,isort)
        jatom=ip(iatom,ig)
        jat=iat_sort_back(jatom)
        ind0=n1pbm*(iat-1)   !! zero index for atom iatom
        jnd0=n1pbm*(jat-1)
        do l=0,lmpb(isort)
          do ie=1,ntle_pb(l,isort)
c ------------- Real components ----------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
              lm0=l+m+1
              t0(lm0)=dreal(f0(jnd0+km))
            enddo
            call rotate1(t0,y0,l,u(1,ig),1)
c ------------- Imaginary components -------------------------------
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
              lm0=l+m+1
              t0(lm0)=dimag(f0(jnd0+km))
            enddo
            call rotate1(t0,y1,l,u(1,ig),1)
            do m=-l,l
              lm=l*(l+1)+m+1
              km=ind_prod(ie,lm,isort)
              lm0=l+m+1
              f(ind0+km)=dcmplx(y0(lm0),y1(lm0))
            enddo
          enddo  !! over ie
        enddo   !! over l
      enddo   !! over iatom
      end
