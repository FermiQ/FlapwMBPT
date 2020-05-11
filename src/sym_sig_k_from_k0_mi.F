      subroutine sym_sig_k_from_k0_mi(k,sig0,sig)
c ------- Finds Sigma at  k = A*k_0 + G_a from Sigma at k_0 ---------      
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k
      complex*16, intent(in) :: sig0(nfun,nbasmpw,npnt)
      complex*16, intent(out) :: sig(nfun,nbasmpw)
      integer :: jbas,jgb,jgb0,jgb1,jrl,jj,ig,j1,gbs_number,k0,
     &           nbas,nb,iatom,i,isort,ia,ib
      real*8 :: v(3),va(3),pi2,tt(3),phase
      complex*16, allocatable :: tmp(:,:),tr(:)
      k0=i_kref(k)
      ig=k_group(k)
      va=gbs(:,g_sym_0(k,ig))
      nbas=nbask(k0)
      nb=nbas/nrel
      sig=(0.d0,0.d0)
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
        do jbas=1,nb
          jgb=indgb(jbas,k)
	      v=gbs(:,jgb)+va  !! G'+G_A
          j1=gbs_number(v)   !! G'+G_A in GBS-list
          j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
          jgb0=iplf_bk(j1,k0)
          if(irel/=2) then
            tmp(i,jbas)=sig0(i,jgb0,k0)
          else if(irel==2) then
            jgb1=jgb0+nb
            do jrl=1,2
              jj=jbas+nb*(jrl-1)
              tmp(i,jj)=tmp(i,jj)+sig0(i,jgb0,k0)*uj(jrl,ig)
     &                           +sig0(i,jgb1,k0)*uj(jrl+2,ig)
            enddo
          endif
        enddo
	  enddo
	  do i=1,nbas
c ------------ Index from the left ------------------------------
        call sym_val_left(ig,tmp(:,i),sig(:,i),1)
	    do iatom=1,natom
	      isort=is(iatom)
	      ia=io_lem(iatom)
	      ib=ia+lfunm(isort)-1
	      sig(ia:ib,i)=sig(ia:ib,i)*tr(iatom)
	    enddo
	  enddo
	  deallocate(tmp,tr)
      end
