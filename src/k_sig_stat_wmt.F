      subroutine k_sig_stat_wmt(it,ind_tau,grs,fif,vloc,s2,ispin,k_pw,n,
     &                          n1,np1,iatom,jatom)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,ind_tau,s2,ispin,n,n1,np1,iatom,jatom
      real*8, intent(in) :: grs(nrel,nfun_red,nfun_red),
     &                      fif(nrel,maxel_red,maxel_red,n_pbmtm_red),
     &                      vloc(n_pbmtm_red,n_pbmtm_red)
      complex*16, intent(out) :: k_pw(n,np1)
      integer :: indi,indj,i,i2,j,l,l2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),gg(:,:)
      allocate(tmp(n,np1))
      indi=io_lem_red(iatom)
      indj=io_lem_red(jatom)
c ------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
      allocate(gg(n,n1))
      if(irel/=2) then
        do j=1,n1
          do i=1,n
            gg(i,j)=grs(1,indi+i-1,indj+j-1)
          enddo
        enddo
      else if(irel==2) then
        do j=1,n1
          do i=1,n
            gg(i,j)=dcmplx(grs(1,indi+i-1,indj+j-1),
     &                     grs(2,indi+i-1,indj+j-1))
          enddo
        enddo
      endif
      tmp=(0.d0,0.d0)
      do i2=1,np1
        do l2=1,n1
          if(irel/=2) cc=fif(1,l2,s2,i2)
          if(irel==2) cc=dcmplx(fif(1,l2,s2,i2),fif(2,l2,s2,i2))
          do l=1,n
            tmp(l,i2)=tmp(l,i2)+gg(l,l2)*cc
          enddo
        enddo
      enddo
      deallocate(gg)
      call dgemm('n','n',2*n,np1,np1,1.d0,tmp,2*n,vloc,n_pbmtm_red,0.d0,
     &           k_pw,2*n)
      deallocate(tmp)
      end
