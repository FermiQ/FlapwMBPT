      subroutine ppm_factors(k,k1,iq,fifim,z1,z2,a1,a2,ppm,gtild,npb,n1,
     &                       n2d,n20,n0)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k1,k,iq,npb,n1,n2d,n20,n0
      real*8, intent(in) :: fifim(maxel,maxel,n_pbmtm,nsort),
     &                      gtild(3,2)
      complex*16, intent(in) :: z1(nfun,n0),z2(nfun,n0),
     &	                      a1(nbasmpw,n0),a2(nbasmpw,n0)
      complex*16, intent(out) :: ppm(n1,n2d,npb)
      integer :: iatom,isort,ind,indp,ndim,ndimp,ia,ib,ic,i,j,km,
     &           ib1,ibas,jbas,id,k0,k10,iq0,id0,km1,ig_tild1,
     &           ig_tild2,igb1,igb2,igb3,ipb,ia0,ib0,ic0,jnd,ia1,ic1,
     &           nbas,nbas1,irl,ist1,ist
      real*8 :: det,det3
      complex*16, allocatable :: tmp(:,:),s(:,:)
      ppm=(0.d0,0.d0)
      k0=i_kref(k)
      k10=i_kref(k1)
      iq0=i_kref(iq)
      allocate(tmp(maxel,n2d))
c ----------------- MT contribution ----------------------------------
      do iatom=1,natom
        isort=is(iatom)
        ind=io_lem(iatom)
        indp=iopb(iatom)-1
        ndim=lfunm(isort)
        ndimp=n_pbmt0(isort)
        do i=1,ndimp
          tmp=(0.d0,0.d0)
          do km1=1,ndim
            j=km1+ind-1
            do km=1,ndim
              do ib=1,n2d
                tmp(km,ib)=tmp(km,ib)+fifim(km,km1,i,isort)*z2(j,n20+ib)
              enddo
            enddo
          enddo
          call zgemm('c','n',n1,n2d,ndim,(1.d0,0.d0),z1(ind,1),nfun,tmp,
     &               maxel,(0.d0,0.d0),ppm(1,1,indp+i),n1)
        enddo   !! over i
      enddo   !! over iatom
      deallocate(tmp)
c ----------------- Interstitial contribution ------------------------
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
      igb1=idnint(det3(gtild(1,1),gbas(1,2),gbas(1,3))/det)
      igb2=idnint(det3(gbas(1,1),gtild(1,1),gbas(1,3))/det)
      igb3=idnint(det3(gbas(1,1),gbas(1,2),gtild(1,1))/det)
      ig_tild1=indplw(igb1,igb2,igb3)
      igb1=idnint(det3(gtild(1,2),gbas(1,2),gbas(1,3))/det)
      igb2=idnint(det3(gbas(1,1),gtild(1,2),gbas(1,3))/det)
      igb3=idnint(det3(gbas(1,1),gbas(1,2),gtild(1,2))/det)
      ig_tild2=indplw(igb1,igb2,igb3)
      allocate(s(nbasmpw,nbasmpw))
      allocate(tmp(nbasmpw,n2d))
      nbas1=nbask(k10)/nrel
      nbas=nbask(k0)/nrel
      do id=1,nplwgw(iq0)
        id0=indpw_gw(id,iq)
        ipb=id+n_pbmt
        ia0=igbs(1,id0)-igbs(1,ig_tild2)+igbs(1,ig_tild1)
        ib0=igbs(2,id0)-igbs(2,ig_tild2)+igbs(2,ig_tild1)
        ic0=igbs(3,id0)-igbs(3,ig_tild2)+igbs(3,ig_tild1)
c ------ We form the overlap matrix -----------------------------
        do jbas=1,nbask(k10)
          jnd=indgb(jbas,k1)
          ia=ia0+igbs(1,jnd)
          ib=ib0+igbs(2,jnd)
          ic=ic0+igbs(3,jnd)
          do ibas=1,nbask(k0)
            ind=indgb(ibas,k)
            ia1=ia-igbs(1,ind)
            ib1=ib-igbs(2,ind)
            ic1=ic-igbs(3,ind)
            j=indplw(ia1,ib1,ic1)
            if(complex_ro) then
              s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
            else
              s(ibas,jbas)=dcmplx(sovr(j),0.d0)
            endif
          enddo
        enddo
        do irl=1,nrel
          ist1=(irl-1)*nbas1
          ist=(irl-1)*nbas
          call zgemm('n','n',nbas,n2d,nbas1,(1.d0,0.d0),s,nbasmpw,
     &               a2(ist1+1,n20+1),nbasmpw,(0.d0,0.d0),tmp,
     &               nbasmpw)
          call zgemm('c','n',n1,n2d,nbas,(1.d0,0.d0),a1(ist+1,1),
     &               nbasmpw,tmp,nbasmpw,(1.d0,0.d0),ppm(1,1,ipb),n1)
        enddo   !! over irl
      enddo   !!! over id
      deallocate(s,tmp)
      end