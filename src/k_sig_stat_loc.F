      subroutine k_sig_stat_loc(it,ind_tau,grs,fif,vloc,s2,k_pw,n,np)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,ind_tau,s2,n,np
      real*8, intent(in) :: grs(nrel,n,n),fif(nrel,n,n,np),
     &                      vloc(np,np)
      complex*16, intent(out) :: k_pw(n,np)
      integer :: i,i2,j,l,l2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),gg(:,:)
      allocate(tmp(n,np))
c ------------------------------------------------------------------
      k_pw=(0.d0,0.d0)
      allocate(gg(n,n))
      if(irel/=2) then
        do j=1,n
          do i=1,n
            gg(i,j)=grs(1,i,j)
          enddo
        enddo
      else if(irel==2) then
        do j=1,n
          do i=1,n
            gg(i,j)=dcmplx(grs(1,i,j),grs(2,i,j))
          enddo
        enddo
      endif
      tmp=(0.d0,0.d0)
      do i2=1,np
        do l2=1,n
          if(irel/=2) cc=fif(1,l2,s2,i2)
          if(irel==2) cc=dcmplx(fif(1,l2,s2,i2),fif(2,l2,s2,i2))
          do l=1,n
            tmp(l,i2)=tmp(l,i2)+gg(l,l2)*cc
          enddo
        enddo
      enddo
      deallocate(gg)
      call dgemm('n','n',2*n,np,np,1.d0,tmp,2*n,vloc,np,0.d0,k_pw,2*n)
      deallocate(tmp)
      end
