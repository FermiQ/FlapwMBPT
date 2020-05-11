	  subroutine sym_sig_red_k_gg(gnew,gold,k0,k,ig,iga)
	  use atom_mod
	  use manager_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
	  implicit none
	  integer, intent(in) :: k,ig,k0,iga
	  complex*16, intent(in) :: gold(nbasmpw_red,nbasmpw_red)
	  complex*16, intent(out) :: gnew(nbasmpw_red,nbasmpw_red)
	  integer :: j,j0,n,ka,n2,ka0,j1,gbs_number,i,i0,i1
	  real*8 :: pi2,v(3),gtild(3),phase
	  complex*16 :: cf,c1,c2,c3,c4
      ka=k_a_from_c(k)
      ka0=k_a_from_c(k0)
      n=nbask_red(ka0)
      n2=n/nrel
	  pi2=pi+pi
      gtild=gbs(:,iga)
	  do j=1,n2
	    j0=indgb_red(j,ka0)  !! G in GBS-list
	    v=gbs(:,j0)+gtild  !! G+G_A
        j1=gbs_number(v)   !! G+G_A in GBS-list
        j1=gbs_sym(j1,ig)   !! A^(-1)*(G+G_A) in GBS-list
        j1=iplf_bk_red(j1,ka)
	    do i=1,n2
	      i0=indgb_red(i,ka0)  !! G in GBS-list
	      v=gbs(:,i0)+gtild  !! G+G_A
          i1=gbs_number(v)   !! G+G_A in GBS-list
          i1=gbs_sym(i1,ig)   !! A^(-1)*(G+G_A) in GBS-list
          i1=iplf_bk_red(i1,ka)
	      v=gbs(:,j0)-gbs(:,i0)  !! G'-G
	      phase=pi2*dot_product(v,shift(:,ig))
	      cf=dcmplx(cos(phase),sin(phase))
	      gnew(i,j)=cf*gold(i1,j1)
          if(irel==2) then
            gnew(n2+i,j)=cf*gold(n2+i1,j1)
            gnew(i,n2+j)=cf*gold(i1,n2+j1)
            gnew(n2+i,n2+j)=cf*gold(n2+i1,n2+j1)
          endif
	    enddo   !! over i
	  enddo  !! over j
      if(irel==2) then
        do i=1,n2
          i1=n2+i
          do j=1,n2
            j1=n2+j
            c1=gnew(i,j)
            c2=gnew(i1,j)
            c3=gnew(i,j1)
            c4=gnew(i1,j1)
            gnew(i,j)=conjg(uj(1,ig))*c1*uj(1,ig)
     &               +conjg(uj(3,ig))*c2*uj(1,ig)
     &               +conjg(uj(1,ig))*c3*uj(3,ig)
     &               +conjg(uj(3,ig))*c4*uj(3,ig)
            gnew(i1,j)=conjg(uj(2,ig))*c1*uj(1,ig)
     &                +conjg(uj(4,ig))*c2*uj(1,ig)
     &                +conjg(uj(2,ig))*c3*uj(3,ig)
     &                +conjg(uj(4,ig))*c4*uj(3,ig)
            gnew(i,j1)=conjg(uj(1,ig))*c1*uj(2,ig)
     &                +conjg(uj(3,ig))*c2*uj(2,ig)
     &                +conjg(uj(1,ig))*c3*uj(4,ig)
     &                +conjg(uj(3,ig))*c4*uj(4,ig)
            gnew(i1,j1)=conjg(uj(2,ig))*c1*uj(2,ig)
     &                 +conjg(uj(4,ig))*c2*uj(2,ig)
     &                 +conjg(uj(2,ig))*c3*uj(4,ig)
     &                 +conjg(uj(4,ig))*c4*uj(4,ig)
          enddo
        enddo
      endif
	  end
