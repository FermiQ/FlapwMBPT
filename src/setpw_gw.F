      subroutine setpw_gw
	  use manager_mod
	  use solid_mod
	  use vertex_mod
      implicit none
      integer :: k,num,idum
      nplw_gw=0
	  do k=1,npnt
	    call g_inside_cutoff(numplw,pnt(1,k),cutoffgw,.false.,.false.,
     &	                     idum,nplwgw(k),nplwpb,nrdiv)
        nplw_gw=max(nplw_gw,nplwgw(k))
      enddo
	  allocate(indpw_gw(nplw_gw,nqdiv))
	  do k=1,npnt
	    call g_inside_cutoff(numplw,pnt(1,k),cutoffgw,.false.,.true.,
     &	                     indpw_gw(1,k),num,nplwpb,nrdiv)
      enddo
      if(mesh_k_c) then
        allocate(nplwgw_red(npnt))
        nplw_gw_red=0
	    do k=1,npnt
	      call g_inside_cutoff(nplwdif,pnt(1,k),cut_pw_red,.false.,
     &	                       .false.,idum,nplwgw_red(k),nplwpb,
     &                         nrdiv_red)
          nplw_gw_red=max(nplw_gw_red,nplwgw_red(k))
        enddo
	    allocate(indpw_gw_red(nplw_gw_red,nqdiv))
	    do k=1,npnt
	      call g_inside_cutoff(nplwdif,pnt(1,k),cut_pw_red,.false.,
     &	                       .true.,indpw_gw_red(1,k),num,nplwpb,
     &                         nrdiv_red)
        enddo
      endif
      end
