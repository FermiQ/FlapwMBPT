      subroutine w_c_from_r_mm(w0,w0r,a,nl,nl1,nat1,nat2,isort,jsort,l,
     &                         l1)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nl,nl1,nat1,nat2,isort,jsort,l,l1
      real*8, intent(in) :: w0r(nl,nl1,nat1,nat2)
      complex*16, intent(in) :: a(2,0:maxpb,natom)
      complex*16, intent(out) :: w0(nl,nl1,nat1,nat2)
      integer :: m,i0,i00,i1,i11,j0,j00,j1,j11,iat,iatom,icase,jat,
     &           jatom,jcase,m1
      do jat=1,nat2
        jatom=iat_sort(jat,jsort)
        j0=ip(jatom,inv_num)
        if(j0==jatom) then
          jcase=1
          j1=j0
        else if(j0<jatom) then
          jcase=2
          j1=jatom
        else if(j0>jatom) then
          jcase=2
          j1=j0
          j0=jatom
        endif
        j00=iat_sort_back(j0)
        j11=iat_sort_back(j1)
        do iat=1,nat1
          iatom=iat_sort(iat,isort)
          i0=ip(iatom,inv_num)
          if(i0==iatom) then
            icase=1
            i1=i0
          else if(i0<iatom) then
            icase=2
            i1=iatom
          else if(i0>iatom) then
            icase=2
            i1=i0
            i0=iatom
          endif
          i00=iat_sort_back(i0)
          i11=iat_sort_back(i1)
          do m1=1,nl1
            do m=1,nl
              if(jcase==1) then
                if(icase==1) then
                  w0(m,m1,i00,j00)=conjg(a(1,l,i0))*a(1,l1,j0)
     &                                             *w0r(m,m1,i00,j00)
                else if(icase==2) then
                  w0(m,m1,i00,j00)=(conjg(a(1,l,i0))*w0r(m,m1,i00,j00)
     &                             +conjg(a(2,l,i0))*w0r(m,m1,i11,j00))
     &                            *a(1,l1,j0)
                  w0(m,m1,i11,j00)=(conjg(a(1,l,i1))*w0r(m,m1,i00,j00)
     &                             +conjg(a(2,l,i1))*w0r(m,m1,i11,j00))
     &                            *a(1,l1,j0)
                endif
              else if(jcase==2) then
                if(icase==1) then
                  w0(m,m1,i00,j00)=conjg(a(1,l,i0))
     &                            *(a(1,l1,j0)*w0r(m,m1,i00,j00)
     &                             +a(2,l1,j0)*w0r(m,m1,i00,j11))
                  w0(m,m1,i00,j11)=conjg(a(1,l,i0))
     &                            *(a(1,l1,j1)*w0r(m,m1,i00,j00)
     &                             +a(2,l1,j1)*w0r(m,m1,i00,j11))
                else if(icase==2) then
                  w0(m,m1,i00,j00)=(conjg(a(1,l,i0))*w0r(m,m1,i00,j00)
     &                             +conjg(a(2,l,i0))*w0r(m,m1,i11,j00))
     &                            *a(1,l1,j0)
     &                            +(conjg(a(1,l,i0))*w0r(m,m1,i00,j11)
     &                             +conjg(a(2,l,i0))*w0r(m,m1,i11,j11))
     &                            *a(2,l1,j0)
                  w0(m,m1,i11,j00)=(conjg(a(1,l,i1))*w0r(m,m1,i00,j00)
     &                             +conjg(a(2,l,i1))*w0r(m,m1,i11,j00))
     &                            *a(1,l1,j0)
     &                            +(conjg(a(1,l,i1))*w0r(m,m1,i00,j11)
     &                             +conjg(a(2,l,i1))*w0r(m,m1,i11,j11))
     &                            *a(2,l1,j0)
                  w0(m,m1,i00,j11)=(conjg(a(1,l,i0))*w0r(m,m1,i00,j00)
     &                             +conjg(a(2,l,i0))*w0r(m,m1,i11,j00))
     &                            *a(1,l1,j1)
     &                            +(conjg(a(1,l,i0))*w0r(m,m1,i00,j11)
     &                             +conjg(a(2,l,i0))*w0r(m,m1,i11,j11))
     &                            *a(2,l1,j1)
                  w0(m,m1,i11,j11)=(conjg(a(1,l,i1))*w0r(m,m1,i00,j00)
     &                             +conjg(a(2,l,i1))*w0r(m,m1,i11,j00))
     &                            *a(1,l1,j1)
     &                            +(conjg(a(1,l,i1))*w0r(m,m1,i00,j11)
     &                             +conjg(a(2,l,i1))*w0r(m,m1,i11,j11))
     &                            *a(2,l1,j1)
                endif
              endif
            enddo  !! over km
          enddo   !! over km1
        enddo   !! over iatom
      enddo   !! over jatom
      end
