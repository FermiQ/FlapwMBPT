      subroutine psi_psi_pw_end(n1,n2,gtild,a1,z1,a2,z2,k,k1,ff,ps3,
     &                          indgbs)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n1,n2,k,k1,indgbs
      real*8, intent(in) :: gtild(3)
      complex*16, intent(in) :: z1(nfun,n1),a1(nbasmpw,n1),z2(nfun,n2),
     &                          a2(nbasmpw,n2),ff(maxel,maxel,natom)
      complex*16, intent(out) :: ps3(n1,n2)
      integer :: iatom,isort,ind,jbas,jnd,ibas,j,ndimb,ind0,gbs_number,
     &           k10,k0
      real*8 :: v(3)
      complex*16, allocatable :: tmp(:,:),s(:,:),tmp1(:,:)
      k0=i_kref(k)
      k10=i_kref(k1)
      ps3=(0.d0,0.d0)
c ----------------- MT contribution --------------------------------
      allocate(tmp(nfun,n2))
      do iatom=1,natom
        isort=is(iatom)
        ndimb=lfunm(isort)
        ind0=io_lem(iatom)
        call zgemm('n','n',ndimb,n2,ndimb,(1.d0,0.d0),ff(1,1,iatom),
     &             maxel,z2(ind0,1),nfun,(0.d0,0.d0),tmp(ind0,1),nfun)
      enddo  !! over iatom
      call zgemm('c','n',n1,n2,nfun,(1.d0,0.d0),z1,nfun,tmp,nfun,
     &           (0.d0,0.d0),ps3,n1)
      deallocate(tmp)
c ----------------- Interstitial contribution ----------------------
c ------ We form the overlap matrix -----------------------------
      allocate(s(nbask(k0),nbask(k10)))
      do jbas=1,nbask(k10)
        jnd=indgb(jbas,k1)
        do ibas=1,nbask(k0)
          ind=indgb(ibas,k)
          v=gbs(:,jnd)-gbs(:,ind)-gtild+gbs(:,indgbs)
          j=gbs_number(v)
          if(complex_ro) then
            s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
          else
            s(ibas,jbas)=dcmplx(sovr(j),0.d0)
          endif
        enddo
      enddo
      allocate(tmp1(nbask(k0),n2))
      call zgemm('n','n',nbask(k0),n2,nbask(k10),(1.d0,0.d0),s,
     &           nbask(k0),a2,nbasmpw,(0.d0,0.d0),tmp1,nbask(k0))
      call zgemm('c','n',n1,n2,nbask(k0),(1.d0,0.d0),a1,nbasmpw,tmp1,
     &           nbask(k0),(1.d0,0.d0),ps3,n1)
      deallocate(s,tmp1)
      ps3=ps3/sqrt(amega)
      end
