      subroutine g_ii_mem_2(gk,g_c,ista,nnn,ind0)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ista,nnn,ind0
      complex*16, intent(in) :: g_c(nbasmpw,nnn,nrel,npnt)
      complex*16, intent(out) :: gk(nbasmpw,nrel,nqdiv)
      integer :: i,jrl,irl,ig,j,iq,iq0,ii,jj,gbs_number,itld,nbas,iba,
     &           krl,jba,j00
      real*8 :: pi2,v(3),phase,phase0
      complex*16 :: cc,yy
      complex*16, allocatable :: xx(:,:)
      pi2=pi+pi
      j=r_pnt_ind(ista)
      allocate(xx(nrel,nrel))
      do iq=1,nqdiv
        ig=k_group(iq)
        iq0=i_kref(iq)
        j00=ind0+ip_r_star(1,ista,ig)
        v=shift(:,ig)+tshift_r(:,j,ig)
        phase0=pi2*dot_product(pnt(:,iq),v)
        itld=g_sym_0(iq,ig)
        nbas=nbask(iq0)/nrel
        do i=1,nbas
          ii=indgb(i,iq)
          v=gbs(:,ii)+gbs(:,itld)
          jj=gbs_number(v)
          jj=gbs_sym(jj,ig)
          jj=iplf_bk(jj,iq0)
          phase=phase0+pi2*dot_product(gbs(:,ii),shift(:,ig))
          cc=dcmplx(cos(phase),-sin(phase))
          if(irel/=2) then
            gk(i,1,iq)=cc*g_c(jj,j00,1,iq0)
          else
            xx=(0.d0,0.d0)
            do irl=1,nrel
              iba=(irl-1)*nbas+jj
              do jrl=1,nrel
                do krl=1,nrel
                  jba=(krl-1)*2+jrl
                  xx(irl,jrl)=xx(irl,jrl)+g_c(iba,j00,krl,iq0)
     &                                   *uj(jba,ig)
                enddo
              enddo
            enddo
            do irl=1,nrel
              iba=(irl-1)*nbas+i
              do jrl=1,nrel
                yy=(0.d0,0.d0)
                do krl=1,nrel
                  jba=(krl-1)*2+irl
                  yy=yy+conjg(uj(jba,ig))*xx(krl,jrl)
                enddo
                gk(iba,jrl,iq)=cc*yy
              enddo
            enddo
          endif
        enddo
      enddo    !! over iq
      deallocate(xx)
      end
