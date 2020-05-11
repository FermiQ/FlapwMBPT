      subroutine p_ii_q_mem_2(pgr,pg,nsta,ista,nnn,ind0)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nsta,ista,nnn,ind0
      complex*16, intent(in) :: pgr(nplw_gw,nqdiv)
      complex*16, intent(inout) :: pg(nplw_gw,nnn,npnt)
      integer :: k,npw,ir1,ir,ig,k0,itld,i,ii,jj,gbs_number
      real*8 :: v(3),phase0,pi2,phase
      complex*16 :: cc
      pi2=pi+pi
      do k=1,npnt
        npw=nplwgw(k)
        do ir1=1,nsta
          ir=ind_r_star(ir1,ista)
          ig=r_group(ir)
          k0=k_sym_0(k,ig)
          v=shift(:,ig)+tshift_r(:,ir,ig)
          phase0=pi2*dot_product(pnt(:,k),v)
          itld=g_sym_0(k,ig)
          do i=1,npw
            ii=indpw_gw(i,k)
            v=gbs(:,ii)+gbs(:,itld)
            jj=gbs_number(v)
            jj=gbs_sym(jj,ig)
            jj=iplf_gk(jj,k0)
            phase=phase0+pi2*dot_product(gbs(:,ii),shift(:,ig))
            cc=dcmplx(cos(phase),-sin(phase))
            pg(i,ind0+ir1,k)=cc*pgr(jj,k0)
          enddo   !! over i
        enddo   !! over ir1
      enddo   !! over k
      end
