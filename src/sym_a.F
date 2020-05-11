      subroutine sym_a(anew,aold,k,ig,pt,n0,neig,nbsk,indb,iplfb)
	  use manager_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: k,ig,n0,neig,nbsk,indb(nbndtm),
     &                       iplfb(nplwbas)
      real*8, intent(in) :: pt(3)
	  complex*16, intent(in) :: aold(n0,neig)
	  complex*16, intent(out) :: anew(n0,neig)
      integer :: j,j0,j1,gbs_number,nb,j2,j3
      real*8 :: v(3),pi2,phase,gtild(3)
	  complex*16 :: cc
	  pi2=pi+pi
      gtild=gbs(:,g_sym_0(k,ig))    !! G_A
      nb=nbsk/nrel
      anew=(0.d0,0.d0)
	  do j=1,nb
	    j0=indb(j)  !! G in GBS-list
	    v=gbs(:,j0)+gtild  !! G+G_A
        j1=gbs_number(v)   !! G+G_A in GBS-list
        j1=gbs_sym(j1,ig)   !! A^(-1)*(G+G_A) in GBS-list
        j1=iplfb(j1)
	    v=pt+gbs(:,j0)  !! k+G
	    phase=-pi2*dot_product(v,shift(:,ig))
	    cc=dcmplx(cos(phase),sin(phase))
        if(irel<2) then
          anew(j,:)=cc*aold(j1,:)
        else if(irel==2) then
          j2=j+nb
          j3=j1+nb
	      anew(j,:)=cc*(conjg(uj(1,ig))*aold(j1,:)
     &                 +conjg(uj(3,ig))*aold(j3,:))
          anew(j2,:)=cc*(conjg(uj(2,ig))*aold(j1,:)
     &                  +conjg(uj(4,ig))*aold(j3,:))
        endif
	  enddo  !! over j
      end
