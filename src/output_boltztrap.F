      subroutine output_boltztrap(ebnd)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: ebnd(nbndf,npnt,nspin)
      integer :: k,info,i,lmax,lmin,ipiv(3),i_len,iatom,isort,ispin,nb,j
      real*8 :: ee,v(3),fc,a(3,3),elim
      elim=0.2d0
      fc=par*0.529177d0
      lmin=10000
      lmax=0
      do ispin=1,nspin
        do k=1,npnt
          nb=n_bnd(k,ispin)
          do i=1,nb
            ee=ebnd(i,k,ispin)-chem_pot
            if(ee>emindos) lmin=min(lmin,i)
            if(ee<emaxdos) lmax=max(lmax,i)
          enddo
        enddo
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'.structure')
          write(3,*)allfile(1:i_len),' basis vectors (Angstroms)'
          do j=1,3
            write(3,'(3(1x,e17.10))')(rbas(i,j)*fc,i=1,3)
          enddo
          write(3,*)natom
          do iatom=1,natom
            isort=is(iatom)
            write(3,'(1x,a4,3(1x,e17.10))')txtel(isort),
     &                                  (tau(i,iatom)*fc,i=1,3)
          enddo
          close(3)
          open(3,file=allfile(1:i_len)//'_DFT.energy')
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_HF.energy')
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_GW.energy')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_scQSGW.energy')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_PSI.energy')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_BSP.energy')
        endif
        write(3,*)' BoltzTrap file, created from FlapwMBPT data'
        write(3,*)npnt,nspin,chem_pot,nelec-(lmin-1)*2/nspin/nrel
        do k=1,npnt
          v=pnt(:,k)
          a=gbas
          call dgesv(3,1,a,3,ipiv,v,3,info)
          write(3,'(3(1x,e17.10),i5)')(v(i),i=1,3),lmax-lmin+1
          do i=lmin,lmax
            write(3,'(2(1x,e17.10))')(ebnd(i,k,ispin),ispin=1,nspin)
          enddo
        enddo
        close(3)
      endif
      end
