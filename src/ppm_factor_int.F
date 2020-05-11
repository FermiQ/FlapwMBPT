      subroutine ppm_factor_int(ipb,k,k1,iq,a1,a2,ppm,gtild,n1,n2)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k1,k,iq,n1,n2,ipb
      real*8, intent(in) :: gtild(3,2)
      complex*16, intent(in) :: a1(nbasmpw,n1),a2(nbasmpw,n2)
      complex*16, intent(out) :: ppm(n1,n2)
      integer :: ind,ibas,jbas,id,k0,k10,iq0,id0,jnd,nbas,nbas1,irl,
     &           ist1,ist,ia,ia0,ia1,ib,ic,ib0,ib1,ic0,ic1,
     &           ig_tild1,ig_tild2,igb1,igb2,igb3,j
      real*8 :: det,det3
      complex*16, allocatable :: tmp(:,:),s(:,:)
      ppm=(0.d0,0.d0)
      k0=i_kref(k)
      k10=i_kref(k1)
      iq0=i_kref(iq)
      nbas1=nbask(k10)/nrel
      nbas=nbask(k0)/nrel
      allocate(s(nbas,nbas1))
      allocate(tmp(nbas,n2))
      id=ipb-n_pbmt
      id0=indpw_gw(id,iq)
      s=(0.d0,0.d0)
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
      igb1=idnint(det3(gtild(1,1),gbas(1,2),gbas(1,3))/det)
      igb2=idnint(det3(gbas(1,1),gtild(1,1),gbas(1,3))/det)
      igb3=idnint(det3(gbas(1,1),gbas(1,2),gtild(1,1))/det)
      ig_tild1=indplw(igb1,igb2,igb3)
      igb1=idnint(det3(gtild(1,2),gbas(1,2),gbas(1,3))/det)
      igb2=idnint(det3(gbas(1,1),gtild(1,2),gbas(1,3))/det)
      igb3=idnint(det3(gbas(1,1),gbas(1,2),gtild(1,2))/det)
      ig_tild2=indplw(igb1,igb2,igb3)
      ia0=igbs(1,id0)-igbs(1,ig_tild2)+igbs(1,ig_tild1)
      ib0=igbs(2,id0)-igbs(2,ig_tild2)+igbs(2,ig_tild1)
      ic0=igbs(3,id0)-igbs(3,ig_tild2)+igbs(3,ig_tild1)
c ------ We form the overlap matrix -----------------------------
      do jbas=1,nbas1
        jnd=indgb(jbas,k1)
        ia=ia0+igbs(1,jnd)
        ib=ib0+igbs(2,jnd)
        ic=ic0+igbs(3,jnd)
        do ibas=1,nbas
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
      s=s/sqrt(amega)
      do irl=1,nrel
        ist1=(irl-1)*nbas1
        ist=(irl-1)*nbas
        call zgemm('n','n',nbas,n2,nbas1,(1.d0,0.d0),s,nbas,
     &             a2(ist1+1,1),nbasmpw,(0.d0,0.d0),tmp,nbas)
        call zgemm('c','n',n1,n2,nbas,(1.d0,0.d0),a1(ist+1,1),
     &             nbasmpw,tmp,nbas,(1.d0,0.d0),ppm,n1)
      enddo   !! over irl
      deallocate(s,tmp)
      end
