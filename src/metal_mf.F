      subroutine metal_mf
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      integer :: ispin,k,n,n1,n2
      integer, allocatable :: n_occ(:,:)
      if(iter>iter_dft+iter_hf) return
      allocate(n_occ(npnt,nspin))
      n_occ=0
      do ispin=1,nspin
        do k=1,npnt
          do i=1,n_bnd(k,ispin)
            x=e_bnd(i,k,ispin)-chem_pot
            if(x<0) n_occ(k,ispin)=i
          enddo
        enddo
      enddo
      metal=.false.
      n1=minval(n_occ(:,1))
      n2=maxval(n_occ(:,1))
      if(n1/=n2) then
        metal=.true.
        goto 1
      endif
      n1=minval(n_occ(:,nspin))
      n2=maxval(n_occ(:,nspin))
      if(n1/=n2) then
        metal=.true.
        goto 1
      endif
      do k=1,npnt
        n=0
        do ispin=1,nspin
          n=n+n_occ(k,ispin)
        enddo
        if(irel/=2.and.nspin==1) n=2*n
        if(n/=nelec) metal=.true.
      enddo
1     deallocate(n_occ)
      if(maswrk) then
        if(metal) write(iun,*)
     &     ' The solid has been identified as METAL'
        if(.not.metal) write(iun,*)
     &     ' The solid has been identified as INSULATOR'
      endif
      end
