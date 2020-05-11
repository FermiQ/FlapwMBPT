      subroutine bravay
	  use atom_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer :: i,iplw,ia,ib,ic
      if(maswrk) then
        write(iun,*)' The vectors of translations :'
        write(iun,'(a3,3f12.7,a3,f12.7,a5)')' A ',(rbas(i,1),i=1,3),
     &    ' = ',sqrt(dot_product(rbas(:,1),rbas(:,1)))*par,' a.u.'
        write(iun,'(a3,3f12.7,a3,f12.7,a5)')' B ',(rbas(i,2),i=1,3),
     &    ' = ',sqrt(dot_product(rbas(:,2),rbas(:,2)))*par,' a.u.'
        write(iun,'(a3,3f12.7,a3,f12.7,a5)')' C ',(rbas(i,3),i=1,3),
     &    ' = ',sqrt(dot_product(rbas(:,3),rbas(:,3)))*par,' a.u.'
      endif
c	---- Get Differences t_iatom - t_jatom  ------------
	  call set_pairs
      call set_lattices
	  call divide_bz
      if(nbndf_bnd_0(2)/=0) call divide_cell_red
c -------------------------------------------------------------
      call getplw
c	---- Get Direct Lattice Vectors encompassing the active polygon --
      call getrbs
c ---------------------------------------------------------------------      
      call neighbors_atoms
c --------------------------------------------------------------
      allocate(indplw(-maxia:maxia,-maxib:maxib,-maxic:maxic))
	  indplw=0
	  do iplw=1,numplw
        ia=igbs(1,iplw)
        ib=igbs(2,iplw)
        ic=igbs(3,iplw)
        indplw(ia,ib,ic)=iplw
      enddo
c --------------------------------------------------------------
	  allocate(i_kref(nqdiv))
	  allocate(k_group(nqdiv))
	  allocate(index_k1(nqdiv))
	  allocate(ip_k(nqdiv,ngroup))
	  allocate(k_star(npnt))
	  allocate(ind_k_star_back(nqdiv))
	  call inf_qnt(npnt,i_kref,k_group,index_k1,pnt,ip_k,k_star,
     &	           k_star_max)
c --------------------------------------------------------------
      if(nbndf_bnd_0(2)/=0) call k_mesh_coarse
c ----------------------------------------------------------------
	  call get_rr
      if(l_double_freq) then
	    call coarse_interp_ac
        call rr_mesh_coarse
        call get_rr_c
      endif
c --------------------------------------------------------------
	  call k_for_band_plot
      call k_sym_intact	
      call getmaxb
      call timel('***** getmaxb finished *************')
      if(maswrk) then
         write(iun,*)' Number of local basis functions ',nlcloc
         write(iun,*)' Max dimension of matrices ',
     &                       nrel*nbndtm+nlcloc
         write(iun,*)' Min dimension of matrices ',
     &                      nrel*nbndtmi+nlcloc
      endif
	  allocate(indgb(nbndtm,nqdiv))
	  allocate(indgb_plot(nbndtm,num_k_plot))
	  call get_indgb
	  call indgb_e
	  allocate(nplwgw(npnt))
	  allocate(iplf_bk(nplwbas,nqdiv))
	  allocate(iplf_gk(nplwpb,nqdiv))
	  allocate(k_sym_0(nqdiv,ngroup))
	  allocate(g_sym_0(nqdiv,ngroup))
	  allocate(gbs_sym(nplwpb,ngroup))
	  call setpw_gw
	  call setpw_gw1
      if(vrt_x_appr(2:2)=='2') call q_long_setup
c --------- For analytical linear K-integration --------------------
      call k_lin_vec
      end
