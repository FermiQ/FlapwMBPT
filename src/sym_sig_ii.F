      subroutine sym_sig_ii(k,k0,ig,s0,s,npw)
c ------- Finds Sigma at  k = A*k_0 + G_a from Sigma at k_0 ---------      
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k,k0,ig,npw
      complex*16, intent(in) :: s0(nrel*npw,nrel*npw)
      complex*16, intent(out) :: s(nrel*npw,nrel*npw)
      integer :: ibas,igb,igb0,jbas,jgb,jgb0,igb1,jgb1,jrl,irl,jj,ii,
     &           j1,gbs_number,nbas
      real*8 :: v(3),va(3)
      complex*16 :: c11,c12,c21,c22,cc
      va=gbs(:,g_sym_0(k,ig))
      nbas=nrel*npw
      s=(0.d0,0.d0)
      do ibas=1,npw
        igb=indpw_gw(ibas,k)
	  v=gbs(:,igb)+va  !! G'+G_A
        j1=gbs_number(v)   !! G'+G_A in GBS-list
        j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
        igb0=iplf_gk(j1,k0)
        do jbas=1,npw
          jgb=indpw_gw(jbas,k)
	    v=gbs(:,jgb)+va  !! G'+G_A
          j1=gbs_number(v)   !! G'+G_A in GBS-list
          j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
          jgb0=iplf_gk(j1,k0)
          if(irel/=2) then
            s(ibas,jbas)=s0(igb0,jgb0)
          else if(irel==2) then
            igb1=igb0+npw
            jgb1=jgb0+npw
            do jrl=1,2
              jj=jbas+npw*(jrl-1)
              do irl=1,2
                ii=ibas+npw*(irl-1)
                c11=uujj(1,1,irl,jrl,ig)
                c21=uujj(2,1,irl,jrl,ig)
                c12=uujj(1,2,irl,jrl,ig)
                c22=uujj(2,2,irl,jrl,ig)
                cc=c11*s0(igb0,jgb0)+c21*s0(igb1,jgb0)+c12*s0(igb0,jgb1)
     &                                                +c22*s0(igb1,jgb1)
                s(ii,jj)=s(ii,jj)+cc
              enddo
            enddo
          endif
        enddo
      enddo   !! over i
      end