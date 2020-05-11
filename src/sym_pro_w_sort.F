      subroutine sym_pro_w_sort(s0,s,ip0,u0,n,npbmt,isort,nat)
      use atom_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ip0(natom),n,npbmt,isort,nat
      real*8, intent(in) :: s0(npbmt,n),u0(maxwig)
      real*8, intent(out) :: s(npbmt,n)
      integer :: iatom,ind0,jatom,jnd0,l,ie,lm,km,l2m,iwig,iat,jat,npb
      s=0.d0
      npb=n_pbmt0(isort)
      do iat=1,nat     !!  over atoms
        iatom=iat_sort(iat,isort)
        ind0=npb*(iat-1)
        jatom=ip0(iatom)
        jat=iat_sort_back(jatom)
        jnd0=npb*(jat-1)
        do l=0,lmpb(isort)
          lm=l*l+1   !! first index for the given L
          l2m=l+l+1
          iwig=l*(2*l-1)*(2*l+1)/3+1
          do ie=1,ntle_pb(l,isort)
            km=ind_prod(ie,lm,isort)
            call dgemm('n','n',l2m,n,l2m,1.d0,u0(iwig),l2m,
     &                 s0(jnd0+km,1),npbmt,0.d0,s(ind0+km,1),npbmt)
          enddo
        enddo   !! over l
      enddo
      end
