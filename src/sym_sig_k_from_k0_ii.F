      subroutine sym_sig_k_from_k0_ii(k,sig0,sig)
c ------- Finds Sigma at  k = A*k_0 + G_a from Sigma at k_0 ---------      
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k
      complex*16, intent(in) :: sig0(nbasmpw,nbasmpw,npnt)
      complex*16, intent(out) :: sig(nbasmpw,nbasmpw)
      integer :: ibas,igb,igb0,jbas,jgb,jgb0,igb1,jgb1,jrl,irl,jj,ii,
     &           ig,j1,gbs_number,k0,nbas,nb
      real*8 :: v(3),va(3)
      complex*16 :: c11,c12,c21,c22,cc
      k0=i_kref(k)
      ig=k_group(k)
      va=gbs(:,g_sym_0(k,ig))
      nbas=nbask(k0)
      nb=nbas/nrel
      sig=(0.d0,0.d0)
      do ibas=1,nb
        igb=indgb(ibas,k)
	    v=gbs(:,igb)+va  !! G'+G_A
        j1=gbs_number(v)   !! G'+G_A in GBS-list
        j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
        igb0=iplf_bk(j1,k0)
        do jbas=1,nb
          jgb=indgb(jbas,k)
	      v=gbs(:,jgb)+va  !! G'+G_A
          j1=gbs_number(v)   !! G'+G_A in GBS-list
          j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
          jgb0=iplf_bk(j1,k0)
          if(irel/=2) then
            sig(ibas,jbas)=sig0(igb0,jgb0,k0)
          else if(irel==2) then
            igb1=igb0+nb
            jgb1=jgb0+nb
            do jrl=1,2
              jj=jbas+nb*(jrl-1)
              do irl=1,2
                ii=ibas+nb*(irl-1)
                c11=uujj(1,1,irl,jrl,ig)
                c21=uujj(2,1,irl,jrl,ig)
                c12=uujj(1,2,irl,jrl,ig)
                c22=uujj(2,2,irl,jrl,ig)
                cc=c11*sig0(igb0,jgb0,k0)
     &            +c21*sig0(igb1,jgb0,k0)
     &            +c12*sig0(igb0,jgb1,k0)
     &            +c22*sig0(igb1,jgb1,k0)
                sig(ii,jj)=sig(ii,jj)+cc
              enddo
            enddo
          endif
        enddo
      enddo   !! over i
      end
