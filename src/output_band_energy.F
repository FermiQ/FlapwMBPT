      subroutine output_band_energy
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      logical :: lkey
      integer :: k,i,nn,i_len,lmin,lmax,n1,n2
      integer, allocatable :: iref(:)
      real*8, allocatable :: d(:)
      lkey=.false.
      if(ubi==' gw') lkey=.true.
      if(ubi==' qp') lkey=.true.
      if(ubi=='psi') lkey=.true.
      if(ubi=='bsp') lkey=.true.
      if(maswrk) then
        i_len=len_trim(allfile)
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_dft.eig')
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_hf.eig')
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_gw_x.eig')
          open(4,file=allfile(1:i_len)//'_gw_qp.eig')
          write(4,*)
     &       ' QP Eigen states relative to the Chemical Potential (eV)'
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_qp_x.eig')
          open(4,file=allfile(1:i_len)//'_qp_qp.eig')
          write(4,*)
     &       ' QP Eigen states relative to the Chemical Potential (eV)'
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_psi_x.eig')
          open(4,file=allfile(1:i_len)//'_psi_qp.eig')
          write(4,*)
     &       ' QP Eigen states relative to the Chemical Potential (eV)'
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_bsp_x.eig')
          open(4,file=allfile(1:i_len)//'_bsp_qp.eig')
          write(4,*)
     &       ' QP Eigen states relative to the Chemical Potential (eV)'
        endif
        write(3,*)
     &       ' Eigen states relative to the Chemical Potential (eV)'
        nn=idnint(nelec*nspin*nrel/2.d0)
        lmin=max(nn-4,1)
        lmax=min(nbndf,lmin+8)
        write(3,'(a8,9i8)')' K\Band ',(i,i=lmin,lmax)
        if(lkey) write(4,'(a8,9i8)')' K\Band ',(i,i=lmin,lmax)
        allocate(d(2*nbndf))
        allocate(iref(2*nbndf))
        do k=1,npnt
          n1=n_bnd(k,1)
          n2=n_bnd(k,nspin)
          if(nspin==1) then
            d(1:n1)=(e_bnd(1:n1,k,1)-chem_pot)*evolt/2.d0
            write(3,'(i5,3x,9f8.4)')k,(d(i),i=lmin,lmax)
            if(lkey) then
              d(1:n1)=(e_qp(1:n1,k,1)-chem_pot)*evolt/2.d0
              write(4,'(i5,3x,9f8.4)')k,(d(i),i=lmin,lmax)
            endif
          else if (nspin==2) then
            d(1:n1)=(e_bnd(1:n1,k,1)-chem_pot)*evolt/2.d0
            d(n1+1:n1+n2)=(e_bnd(1:n2,k,2)-chem_pot)*evolt/2.d0
            call ordering_all(n1+n2,d,iref)
            write(3,'(i5,3x,9f8.4)')k,(d(iref(i)),i=lmin,lmax)
            if(lkey) then
              d(1:n1)=(e_qp(1:n1,k,1)-chem_pot)*evolt/2.d0
              d(n1+1:n1+n2)=(e_qp(1:n2,k,2)-chem_pot)*evolt/2.d0
              call ordering_all(n1+n2,d,iref)
              write(4,'(i5,3x,9f8.4)')k,(d(iref(i)),i=lmin,lmax)
            endif
          endif
        enddo   !!! over k
        deallocate(d,iref)
        close(3)
        if(lkey) close(4)
      endif
      end
