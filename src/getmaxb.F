      subroutine getmaxb
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer :: ipnt,nbndt,ierr,k,idum
	  real*8 :: cut
	  allocate(nbask(npnt),STAT=ierr)
	  allocate(nbask_plot(num_k_plot),STAT=ierr)
	  allocate(n_bnd(npnt,nspin_0),STAT=ierr)
      nbndtm=0
      nbndtmi=1000000
      do ipnt=1,npnt
        call g_inside_cutoff(nplwbas,pnt(1,ipnt),cutoff,.false.,.false.,
     &                       idum,nbndt,nplwbas,nrdiv)
	    if(nbndf==0) n_bnd(ipnt,:)=nrel*nbndt+nlcloc
	    nbask(ipnt)=nrel*nbndt
        nbndtm=max(nbndtm,nbndt)
        nbndtmi=min(nbndtmi,nbndt)
      enddo
      do k=1,num_k_plot
        call g_inside_cutoff(nplwbas,k_plot(1,k),cutoff,.false.,.false.,
     &                       idum,nbndt,nplwbas,nrdiv)
	    nbask_plot(k)=nrel*nbndt
        nbndtm=max(nbndtm,nbndt)
      enddo
      nbasmpw=nrel*nbndtm
	  nbasm=nbasmpw+nlcloc
	  if(nbndf==0) nbndf=nbasm
	  if(nbndf>nbasm) nbndf=nbasm
	  allocate(e_bnd(nbndf,npnt,nspin_0))
c --------------------------------------------------------------------	
      if(mesh_k_c) then
	    allocate(nbask_red(npnt),STAT=ierr)
	    cut=min(cutoff,cut_bnd_red)
        nbasmpw_red=0
        do ipnt=1,npnt
	      call g_inside_cutoff(nplwbas,pnt(1,ipnt),cut,.false.,
     &                         .false.,idum,nbndt,nplwbas,nrdiv_red)
	      nbask_red(ipnt)=nrel*nbndt
          nbasmpw_red=max(nbasmpw_red,nbndt)
        enddo
        nbasmpw_red=nrel*nbasmpw_red
        if(maswrk)write(iun,'(a51,i4,a1,i4)')
     &     ' Min/Max number of (sG)-waves for low energy bands:',
     &              minval(nbask_red),'/',maxval(nbask_red)
      endif
      end
