      subroutine get_rr_c
      use atom_mod
      use manager_mod
      use parallel_mod
      use sitegen_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer :: n1,n2,n3,ind,ind0,i1,i2,i3,ig,iatom,jatom,iatom0,
     &           jatom0
      real*8 :: v(3),v1(3),v2(3)
      integer, allocatable :: iflag(:,:,:)
      real*8, allocatable :: rr_pntc(:,:)
      allocate(rr_pntc(3,nqdiv_c))
      allocate(iflag(nqdiv_c,natom,natom))
c  set divisions along x,y,z
      allocate(rr_c3(3,nqdiv_c))
      n1=ndiv_c(1)-1
      n2=ndiv_c(2)-1
      n3=ndiv_c(3)-1
      ind0=0
      do i3=0,n3
        do i2=0,n2
          do i1=0,n1
            ind0=ind0+1
            rr_pntc(1,ind0)=i1
            rr_pntc(2,ind0)=i2
            rr_pntc(3,ind0)=i3
            rr_c3(:,ind0)=i1*rbas(:,1)+i2*rbas(:,2)+i3*rbas(:,3)
          enddo
        enddo
      enddo
      allocate(r0_pair_c(nqdiv_c,natom,natom))
      allocate(ig_pair_c(nqdiv_c,natom,natom))
      allocate(nrr_red_c(3,nqdiv_c*natom*natom))
      nrr_max_c=0
      nrr_red_c=0
      r0_pair_c=0
      iflag=0
      do jatom0=1,natom
        do iatom0=1,natom
          do ind0=1,nqdiv_c
            if(iflag(ind0,iatom0,jatom0)/=0) cycle
            nrr_max_c=nrr_max_c+1
            nrr_red_c(1,nrr_max_c)=ind0
            nrr_red_c(2,nrr_max_c)=iatom0
            nrr_red_c(3,nrr_max_c)=jatom0
            i1=rr_pntc(1,ind0)
            i2=rr_pntc(2,ind0)
            i3=rr_pntc(3,ind0)
            v=i1*rbas(:,1)+i2*rbas(:,2)+i3*rbas(:,3)
c  rotate given point
            do ig=1,ngroup
              do i1=1,natom
                if(ip(i1,ig)==iatom0) then
                  iatom=i1
                  exit
                endif
              enddo
              do i1=1,natom
                if(ip(i1,ig)==jatom0) then
                  jatom=i1
                  exit
                endif
              enddo
              v1=v+tau(:,iatom0)-tau(:,jatom0)
              call rotate(v1(1),v1(2),v1(3),v2(1),v2(2),v2(3),u(2,ig),2)
              v2=v2-tau(:,iatom)+tau(:,jatom)
              call zone1_number(v2,gbas,ndiv_c,ind)
              if(iflag(ind,iatom,jatom)==1) cycle
              iflag(ind,iatom,jatom)=1
              r0_pair_c(ind,iatom,jatom)=nrr_max_c
              ig_pair_c(ind,iatom,jatom)=ig
            enddo
          enddo
        enddo  !! over iatom0
      enddo  !! over jatom0
      deallocate(iflag,rr_pntc)
      if(maswrk) write(iun,*)' NRR_MAX_C = ',nrr_max_c
      end
