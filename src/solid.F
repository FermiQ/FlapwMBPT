	  subroutine solid
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
      implicit none
	  integer :: maxbv,maxbj,maxw2,ind
	  ind_tens33(1,1)=1
      ind_tens33(2,1)=4
      ind_tens33(3,1)=5
      ind_tens33(1,2)=4
      ind_tens33(2,2)=2
      ind_tens33(3,2)=6
      ind_tens33(1,3)=5
      ind_tens33(2,3)=6
      ind_tens33(3,3)=3
      maxbv=max(maxb,maxpb,1)
      maxbj=2*maxbv+1
      maxw2=max(maxbv,1)
      maxwig=(maxbv+1)*(2*maxbv+1)*(2*maxbv+3)/3
      maxwig_j=4*(maxbv+1)*(2*(maxbv+1)**2+1)/3
      allocate(ifia(148))
      allocate(shift(3,148))
      allocate(u(maxwig,148))
      allocate(u_opt(3,3,148))
      allocate(ag(3,3,148))
      allocate(agm(3,3,148))
      allocate(uj(maxwig_j,148))
      call get_bz
      call group(maxbv,maxbj,maxw2)
      call geometry
	  call set_solid
      call memory_solid(maxbv,maxbj,maxw2)
! ------- Real mesh checkings ------------------------------------
      call r_mesh_checking(mdiv)
	  if(maswrk) write(iun,'(a35,3i5)')
     &  ' R_MESH_CHECKING: MDIV is set to be',(mdiv(ind),ind=1,3)
      call r_mesh_checking(nrdiv)
	  if(maswrk) write(iun,'(a36,3i5)')
     &  ' R_MESH_CHECKING: NRDIV is set to be',(nrdiv(ind),ind=1,3)
      call r_mesh_checking(nrdiv_red)
	  if(maswrk) write(iun,'(a40,3i5)')
     &  ' R_MESH_CHECKING: NRDIV_RED is set to be',
     &                                     (nrdiv_red(ind),ind=1,3)
! -----------------------------------------------------------------     
      call bravay
c      call transport_boltz
      call symlmm
      if(maxnsymb.eq.0) maxnsymb=1
	  call getmts
	  call radmsh
      maxnrad=maxval(nrad)
	  call start_par_solid
      call getind
c -------------------------------------------------------------------
	  call memory1_solid
	  call freeatom
      call timel('**** Freeatom finished *************')
	  call neighbor_atm
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      allocate(sovr(numplw))
      if(complex_ro) allocate(sovi(numplw))
      call warp
c ----------------------------------------------------------------
	  call memory_flapw
	  call solid_scf
	  call memory_close
      end
