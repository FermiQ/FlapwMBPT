      subroutine sa_prep(sa0)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(out) :: sa0(nplw_gw,nbndf,ndim3_k(me_k+1),nrel,
     &                               nspin)
      integer :: ispin,ind_k,k,ib,i,i0,ia,ibas0,ic,ind,irl,j,j0,nbas
      complex*16, allocatable :: s(:,:)
      allocate(s(nplw_gw,nbndtm))
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        nbas=nbask(k)/nrel
c ------ We form the overlap matrix -----------------------------
        do j=1,nbas
          j0=indgb(j,k)
          do i=1,nplwgw(k)
            i0=indpw_gw(i,k)
            ia=igbs(1,j0)-igbs(1,i0)
            ib=igbs(2,j0)-igbs(2,i0)
            ic=igbs(3,j0)-igbs(3,i0)
            ind=indplw(ia,ib,ic)
            if(complex_ro) then
              s(i,j)=dcmplx(sovr(ind),sovi(ind))
            else
              s(i,j)=dcmplx(sovr(ind),0.d0)
            endif
          enddo
        enddo
        do ispin=1,nspin
c -------- Forming S*A ---------------------------------------------
          do irl=1,nrel
            ibas0=(irl-1)*nbas
            call zgemm('n','n',nplwgw(k),n_bnd(k,ispin),nbas,
     &                 (1.d0,0.d0),s,nplw_gw,
     &                 ev_bnd(ibas0+1,1,ind_k,ispin),nbasmpw,
     &                 (0.d0,0.d0),sa0(1,1,ind_k,irl,ispin),nplw_gw)
          enddo
        enddo   !! over ispin
      enddo   !! over ind_k
      deallocate(s)
      end
