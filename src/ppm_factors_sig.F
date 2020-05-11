      subroutine ppm_factors_sig(k,k1,iq,fif_full,z_full,a_full,ppm,
     &                           gtild,num_bands)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k1,k,iq,num_bands(npnt)
      real*8, intent(in) :: fif_full(nrel,maxel**2,n_pbmtm_red,nsort),
     &                      gtild(3)
      complex*16, intent(in) :: z_full(nfun,nbndf_bnd,nqdiv_c),
     &                        a_full(nbasmpw,nbndf_bnd,nqdiv_c)
      complex*16, intent(out) :: ppm(n_pbtot_red,nbndf_bnd,nbndf_bnd)
      integer :: iatom,isort,ind,indp,ndim,ndimp,ia,ib,ic,i,j,km,ij,
     &           ib1,ibas,jbas,id,k0,k10,iq0,n0,n1,id0,km1,ig_tild,
     &           igb1,igb2,igb3,ipb,ia0,ib0,ic0,jnd,ia1,ic1,ii,
     &           kf0,kf10,kf1,kf,iqa,iqa0,nbas,nbas1,irl,ist1,ist
      real*8 :: det,det3
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),s(:,:),tt(:,:,:)
      allocate(tt(nbndf_bnd,nbndf_bnd,n_pbtot_red))
      tt=(0.d0,0.d0)
      ppm=(0.d0,0.d0)
      k0=i_kref_c(k)
      k10=i_kref_c(k1)
      kf=k_a_from_c(k)
      kf1=k_a_from_c(k1)
      kf0=k_a_from_c(k0)
      kf10=k_a_from_c(k10)
      iq0=i_kref_c(iq)
      iqa0=k_a_from_c(iq0)
      iqa=k_a_from_c(iq)
      n0=num_bands(kf0)
      n1=num_bands(kf10)
      allocate(tmp(nbndf_bnd,maxel))
c ----------------- MT contribution ----------------------------------
      do iatom=1,natom
        isort=is(iatom)
        ind=io_lem(iatom)
        indp=iopb_red(iatom)-1
        ndim=lfunm(isort)
        ndimp=n_pbmt0_red(isort)
        do i=1,ndimp
          tmp=(0.d0,0.d0)
          ij=0
          do km1=1,ndim
            do km=1,ndim
              ij=ij+1
              if(irel/=2) cc=fif_full(1,ij,i,isort)
              if(irel==2) cc=dcmplx(fif_full(1,ij,i,isort),
     &                              fif_full(2,ij,i,isort))
              do ib=1,n0
                tmp(ib,km1)=tmp(ib,km1)+cc*conjg(z_full(km+ind-1,ib,k))
              enddo
            enddo
          enddo
          call zgemm('n','n',n0,n1,ndim,(1.d0,0.d0),tmp,nbndf_bnd,
     &             z_full(ind,1,k1),nfun,(0.d0,0.d0),
     &               tt(1,1,indp+i),nbndf_bnd)
        enddo   !! over i
      enddo   !! over iatom
      deallocate(tmp)
c ----------------- Interstitial contribution ------------------------
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
      igb1=idnint(det3(gtild,gbas(1,2),gbas(1,3))/det)
      igb2=idnint(det3(gbas(1,1),gtild,gbas(1,3))/det)
      igb3=idnint(det3(gbas(1,1),gbas(1,2),gtild)/det)
      ig_tild=indplw(igb1,igb2,igb3)
      allocate(s(nbasmpw,nbasmpw))
      allocate(tmp(nbasmpw,nbndf_bnd))
      nbas1=nbask(kf10)/nrel
      nbas=nbask(kf0)/nrel
      do id=1,nplwgw_red(iqa0)
        id0=indpw_gw_red(id,iqa)
        ipb=id+n_pbmt_red
        ia0=igbs(1,id0)-igbs(1,ig_tild)
        ib0=igbs(2,id0)-igbs(2,ig_tild)
        ic0=igbs(3,id0)-igbs(3,ig_tild)
c ------ We form the overlap matrix -----------------------------
        do jbas=1,nbas1
          jnd=indgb(jbas,kf1)
          ia=ia0+igbs(1,jnd)
          ib=ib0+igbs(2,jnd)
          ic=ic0+igbs(3,jnd)
          do ibas=1,nbas
            ind=indgb(ibas,kf)
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
          call zgemm('n','n',nbas,n1,nbas1,(1.d0,0.d0),s,nbasmpw,
     &               a_full(ist1+1,1,k1),nbasmpw,(0.d0,0.d0),tmp,
     &               nbasmpw)
          call zgemm('c','n',n0,n1,nbas,(1.d0,0.d0),a_full(ist+1,1,k),
     &               nbasmpw,tmp,nbasmpw,(1.d0,0.d0),tt(1,1,ipb),
     &               nbndf_bnd)
        enddo   !! over irl
      enddo   !!! over id
      do ii=1,n_pbtot_red
        do j=1,nbndf_bnd
          ppm(ii,j,:)=tt(:,j,ii)
        enddo
      enddo
      deallocate(s,tmp,tt)
      end
