      subroutine get_rr
	  use atom_mod
	  use manager_mod
	  use parallel_mod
      use sitegen_mod
	  use solid_mod
      implicit none
      integer :: n1,n2,n3,ind,ind0,i1,i2,i3,ig,iatom,jatom,iatom0,
     &           jatom0
      real*8 :: v(3),v1(3),v2(3)
	  integer, allocatable :: iflag(:,:,:)
      allocate(rr_pnt(3,nqdiv))
      allocate(rr_3(3,nqdiv))
	  allocate(iflag(nqdiv,natom,natom))
c  set divisions along x,y,z
      n1=ndiv(1)-1
      n2=ndiv(2)-1
      n3=ndiv(3)-1
	  ind0=0
      do i3=0,n3
        do i2=0,n2
          do i1=0,n1
	        ind0=ind0+1
	        rr_pnt(1,ind0)=i1
	        rr_pnt(2,ind0)=i2
	        rr_pnt(3,ind0)=i3
            rr_3(:,ind0)=i1*rbas(:,1)+i2*rbas(:,2)+i3*rbas(:,3)
          enddo
        enddo
      enddo
      allocate(r0_pair(nqdiv,natom,natom))
      allocate(ig_pair(nqdiv,natom,natom))
      allocate(nrr_red(3,nqdiv*natom*natom))
	  nrr_max=0
	  nrr_red=0
	  r0_pair=0
	  iflag=0
	  do jatom0=1,natom
	    do iatom0=1,natom
	      do ind0=1,nqdiv
	        if(iflag(ind0,iatom0,jatom0)/=0) cycle
	        nrr_max=nrr_max+1
	        nrr_red(1,nrr_max)=ind0
	        nrr_red(2,nrr_max)=iatom0
	        nrr_red(3,nrr_max)=jatom0
	        i1=rr_pnt(1,ind0)
	        i2=rr_pnt(2,ind0)
	        i3=rr_pnt(3,ind0)
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
	          call zone1_number(v2,gbas,ndiv,ind)
	          if(iflag(ind,iatom,jatom)==1) cycle
	          iflag(ind,iatom,jatom)=1
	          r0_pair(ind,iatom,jatom)=nrr_max
	          ig_pair(ind,iatom,jatom)=ig
	        enddo
	      enddo
	    enddo  !! over iatom0
	  enddo  !! over jatom0
	  deallocate(iflag)
	  if(maswrk) write(iun,*)' NRR_MAX = ',nrr_max
      end
