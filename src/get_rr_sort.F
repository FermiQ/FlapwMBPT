      subroutine get_rr_sort(isort,jsort,nat1,nat2,r0_pairs,ig_pairs,
     &                       nrr_reds,nrr_maxs)
      use atom_mod
      use manager_mod
      use parallel_mod
      use sitegen_mod
      use solid_mod
      implicit none
      integer, intent(in) :: isort,jsort,nat1,nat2
      integer, intent(out) :: r0_pairs(nqdiv,nat1,nat2),nrr_maxs,
     &                        ig_pairs(nqdiv,nat1,nat2),
     &                        nrr_reds(3,nqdiv*nat1*nat2)
      integer :: ind,ind0,i1,i2,i3,ig,iatom,jatom,iatom0,jatom0,i11,
     &           jat0,iat0,iat,iatom1,jat,jatom1
      real*8 :: v(3),v1(3),v2(3)
      integer, allocatable :: iflag(:,:,:)
      allocate(iflag(nqdiv,nat1,nat2))
      nrr_maxs=0
      nrr_reds=0
      r0_pairs=0
      iflag=0
      do jatom0=1,nat2
        jat0=iat_sort(jatom0,jsort)
        do iatom0=1,nat1
          iat0=iat_sort(iatom0,isort)
          do ind0=1,nqdiv
            if(iflag(ind0,iatom0,jatom0)/=0) cycle
            nrr_maxs=nrr_maxs+1
            nrr_reds(1,nrr_maxs)=ind0
            nrr_reds(2,nrr_maxs)=iatom0
            nrr_reds(3,nrr_maxs)=jatom0
            i1=rr_pnt(1,ind0)
            i2=rr_pnt(2,ind0)
            i3=rr_pnt(3,ind0)
            v=i1*rbas(:,1)+i2*rbas(:,2)+i3*rbas(:,3)
c  rotate given point
            do ig=1,ngroup
              do i1=1,nat1
                iatom1=iat_sort(i1,isort)
                iatom1=ip(iatom1,ig)
                i11=iat_sort_back(iatom1)
                if(i11==iatom0) then
                  iatom=i1
                  iat=iat_sort(iatom,isort)
                  exit
                endif
              enddo
              do i1=1,nat2
                jatom1=iat_sort(i1,jsort)
                jatom1=ip(jatom1,ig)
                i11=iat_sort_back(jatom1)
                if(i11==jatom0) then
                  jatom=i1
                  jat=iat_sort(jatom,jsort)
                  exit
                endif
              enddo
              v1=v+tau(:,iat0)-tau(:,jat0)
              call rotate(v1(1),v1(2),v1(3),v2(1),v2(2),v2(3),u(2,ig),2)
              v2=v2-tau(:,iat)+tau(:,jat)
              call zone1_number(v2,gbas,ndiv,ind)
              if(iflag(ind,iatom,jatom)==1) cycle
              iflag(ind,iatom,jatom)=1
              r0_pairs(ind,iatom,jatom)=nrr_maxs
              ig_pairs(ind,iatom,jatom)=ig
            enddo
          enddo
        enddo  !! over iatom0
      enddo  !! over jatom0
      deallocate(iflag)
      end
