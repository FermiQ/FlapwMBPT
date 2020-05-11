      subroutine sym_sig_mi(k,k0,s0,s,ig,npw)
c ------- Finds Sigma at  k = A*k_0 + G_a from Sigma at k_0 ---------      
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k,ig,npw
      complex*16, intent(in) :: s0(nfun,nrel*npw)
      complex*16, intent(out) :: s(nfun,nrel*npw)
      integer :: jbas,jgb,jgb0,jgb1,jrl,jj,j1,gbs_number,k0,nbas,iatom,
     &           i,isort,ia,ib
      real*8 :: v(3),va(3),pi2,tt(3),phase
      complex*16, allocatable :: tmp(:,:),tr(:)
      va=gbs(:,g_sym_0(k,ig))
      nbas=nrel*npw
	pi2=pi+pi
      allocate(tmp(nfun,nbas))
      allocate(tr(natom))
	do iatom=1,natom
	  tt=tshift(:,iatom,ig)
	  phase=pi2*dot_product(pnt(:,k),tt)
	  tr(iatom)=dcmplx(cos(phase),sin(phase))
	enddo
	tmp=(0.d0,0.d0)
	do i=1,nfun
c ------------ Index from the right  ------------------------------
        do jbas=1,npw
          jgb=indpw_gw(jbas,k)
	    v=gbs(:,jgb)+va  !! G'+G_A
          j1=gbs_number(v)   !! G'+G_A in GBS-list
          j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
          jgb0=iplf_gk(j1,k0)
          if(irel/=2) then
            tmp(i,jbas)=s0(i,jgb0)
          else if(irel==2) then
            jgb1=jgb0+npw
            do jrl=1,2
              jj=jbas+npw*(jrl-1)
              tmp(i,jj)=tmp(i,jj)+s0(i,jgb0)*uj(jrl,ig)
     &                           +s0(i,jgb1)*uj(jrl+2,ig)
            enddo
          endif
        enddo
      enddo
      s=(0.d0,0.d0)
	do i=1,nbas
c ------------ Index from the left ------------------------------
        call sym_val_left(ig,tmp(:,i),s(:,i),1)
	  do iatom=1,natom
	    isort=is(iatom)
	    ia=io_lem(iatom)
	    ib=ia+lfunm(isort)-1
	    s(ia:ib,i)=s(ia:ib,i)*tr(iatom)
	  enddo
	enddo
	deallocate(tmp,tr)
      end