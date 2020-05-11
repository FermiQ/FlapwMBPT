      subroutine ppm_factors0(k,k1,iq,fifim,z1,z2,a1,a2,ppm,gtild,npb,
     &                        n1,n2d,n20,n0)
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
     &                          a1(nbasmpw,n0),a2(nbasmpw,n0)
      complex*16, intent(out) :: ppm(n1,n2d,npb)
      integer :: iatom,isort,ind,indp,ndim,ndimp,ib,i,j,km,ibas,jbas,id,
     &           k0,k10,iq0,id0,km1,ipb,jnd,nbas,nbas1,irl,ist1,ist,
     &           gbs_number
      real*8 :: v(3),rc
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
c ----------------- Interstitial contribution ----------------------
      nbas1=nbask(k10)/nrel
      nbas=nbask(k0)/nrel
      allocate(s(nbas,nbas1))
      allocate(tmp(nbas,n2d))
      rc=1.d0/sqrt(amega)
      do id=1,nplwgw(iq0)
        id0=indpw_gw(id,iq)
        ipb=id+n_pbmt
        s=(0.d0,0.d0)
        do jbas=1,nbas1
          jnd=indgb(jbas,k1)
          do ibas=1,nbas
            ind=indgb(ibas,k)
            v=gbs(:,ind)+gtild(:,2)-gbs(:,jnd)
            if(gbs_number(v)==id0) s(ibas,jbas)=rc
          enddo
        enddo
        do irl=1,nrel
          ist1=(irl-1)*nbas1
          ist=(irl-1)*nbas
          call zgemm('n','n',nbas,n2d,nbas1,(1.d0,0.d0),s,nbas,
     &               a2(ist1+1,n20+1),nbasmpw,(0.d0,0.d0),tmp,nbas)
          call zgemm('c','n',n1,n2d,nbas,(1.d0,0.d0),a1(ist+1,1),
     &               nbasmpw,tmp,nbas,(1.d0,0.d0),ppm(1,1,ipb),n1)
        enddo   !! over irl
      enddo   !!! over id
      deallocate(s,tmp)
      end
