	  subroutine sym_sig_red_k_gm(gnew,gold,k0,k,ig,iga)
	  use atom_mod
	  use manager_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
	  implicit none
	  integer, intent(in) :: k,ig,k0,iga
	  complex*16, intent(in) :: gold(nbasmpw_red,nfun_red)
	  complex*16, intent(out) :: gnew(nbasmpw_red,nfun_red)
	  integer :: iatom,isort,ind0,ind1,j,j0,n,n2,ka0,j1,gbs_number,i,
     &           ka
	  real*8 :: pi2,fas,v(3),gtild(3),phase
	  complex*16 :: cf,c1,c2
	  complex*16, allocatable :: tmp(:,:)
      ka=k_a_from_c(k)
      ka0=k_a_from_c(k0)
      n=nbask_red(ka0)
      n2=n/nrel
	  pi2=pi+pi
	  allocate(tmp(n,nfun_red))
	  do j=1,n
	    call sym_left_red(ig,gold(j,:),tmp(j,:))
      enddo
      gtild=gbs(:,iga)
	  do j=1,n2
	    j0=indgb_red(j,ka0)  !! G in GBS-list
	    v=gbs(:,j0)+gtild  !! G+G_A
        j1=gbs_number(v)   !! G+G_A in GBS-list
        j1=gbs_sym(j1,ig)   !! A^(-1)*(G+G_A) in GBS-list
        j1=iplf_bk_red(j1,ka)
	    v=pnt_c(:,k0)+gbs(:,j0)  !! k+G
	    phase=-pi2*dot_product(v,shift(:,ig))
	    do iatom=1,natom
	      isort=is(iatom)
	      ind0=io_lem_red(iatom)
	      ind1=ind0+lfunm_red(isort)-1
	      v=tshift(:,iatom,ig)
 	      fas=phase-pi2*dot_product(pnt_c(:,k0),v)
	      cf=dcmplx(cos(fas),sin(fas))
	      gnew(j,ind0:ind1)=cf*tmp(j1,ind0:ind1)
          if(irel==2) gnew(n2+j,ind0:ind1)=cf*tmp(n2+j1,ind0:ind1) 
	    enddo   !! over iatom_c
	  enddo  !! over j
	  deallocate(tmp)
      if(irel==2) then
        do i=1,nfun_red
          do j=1,n2
            j1=n2+j
            c1=gnew(j,i)
            c2=gnew(j1,i)
            gnew(j,i)=conjg(uj(1,ig))*c1+conjg(uj(3,ig))*c2
            gnew(j1,i)=conjg(uj(2,ig))*c1+conjg(uj(4,ig))*c2
          enddo
        enddo
      endif
	  end
