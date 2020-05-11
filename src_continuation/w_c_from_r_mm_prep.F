      subroutine w_c_from_r_mm_prep(a)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(out) :: a(2,0:maxpb,natom,ndim3_k(me_k+1))
      integer :: npw,jcase,icase,jatom,j0,j1,i0,i1,iatom,isort,jsort,l,
     &           l1,jnd0,jnd1,ind0,ind1,ii1,jj1,lm1,lm,nn1,j,i,ii2,jj2,
     &           km,km1,lget,ind_k,k
      real*8 :: phase,pi2
      complex*16 :: c1,c3
      pi2=pi+pi
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        do iatom=1,natom
          isort=is(iatom)
          i0=ip(iatom,inv_num)
          if(i0==iatom) then
            phase=-pi2*dot_product(pnt(:,k),tau(:,iatom))
            c1=dcmplx(cos(phase),sin(phase))
            do l=0,lmpb(isort)
              a(1,l,iatom,ind_k)=(0.d0,-1.d0)**l*c1
            enddo
          else if(i0<iatom) then
            phase=-pi2*dot_product(pnt(:,k),tau(:,iatom))
            c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
            do l=0,lmpb(isort)
              c3=(0.d0,-1.d0)**l
              a(1,l,iatom,ind_k)=c3*c1
              a(2,l,iatom,ind_k)=c3*(0.d0,-1.d0)*c1
            enddo
          else if(i0>iatom) then
            phase=-pi2*dot_product(pnt(:,k),tau(:,iatom))
            c1=dcmplx(cos(phase),sin(phase))/sqrt(2.d0)
            do l=0,lmpb(isort)
              c3=(0.d0,-1.d0)**l
              a(1,l,iatom,ind_k)=c3*c1
              a(2,l,iatom,ind_k)=-c3*(0.d0,-1.d0)*c1
            enddo
          endif
        enddo
      enddo
      end
