      subroutine sig_1_stat_loc(it,ind_tau,t_pw,v00,fif,s2,grs,sig_stat,
     &                          n,np)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,ind_tau,s2,n,np
      real*8, intent(in) :: grs(nrel,n,n,0:n_tau),
     &                      fif(nrel,n,n,np),
     &                      v00(np,np)
      complex*16, intent(inout) :: sig_stat(n)
      complex*16, intent(out) :: t_pw(np,n)
      integer :: i_tau,i1_tau,i,ii,j
      complex*16 :: cc
      complex*16, allocatable :: k_pw(:,:),w(:,:),tmp(:,:),tmp1(:,:),
     &                           gg(:,:)
      i_tau=me_t*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau
      i1_tau=n_tau-i_tau
      allocate(k_pw(n,np))
      call k_sig_stat_loc(it,ind_tau,grs(1,1,1,i_tau),fif,v00,s2,k_pw,n,
     &                    np)
      if(irel/=2) then
        call tsig_loc(grs(1,1,1,i1_tau),k_pw,t_pw,fif,n,np)
      else
        call tsig_loc_c(grs(1,1,1,i1_tau),k_pw,t_pw,fif,n,np)
      endif
      deallocate(k_pw)
      allocate(tmp(np,n))
      allocate(tmp1(np,n))
      allocate(w(np,np))
      w=v00
      call zgemm('n','n',np,n,np,(1.d0,0.d0),w,np,t_pw,np,(0.d0,0.d0),
     &           tmp,np)
      deallocate(w)
      allocate(gg(n,n))
      if(irel/=2) then
        do j=1,n
          do i=1,n
            gg(i,j)=grs(1,i,j,i_tau)
          enddo
        enddo
      else if(irel==2) then
        do j=1,n
          do i=1,n
            gg(i,j)=dcmplx(grs(1,i,j,i_tau),
     &                     grs(2,i,j,i_tau))
          enddo
        enddo
      endif
      call zgemm('n','t',np,n,n,(1.d0,0.d0),tmp,np,gg,n,(0.d0,0.d0),
     &           tmp1,np)
      deallocate(gg,tmp)
      do ii=1,n
        do j=1,n
          do i=1,np
            if(irel/=2) cc=fif(1,ii,j,i)
            if(irel==2) cc=dcmplx(fif(1,ii,j,i),fif(2,ii,j,i))
            sig_stat(ii)=sig_stat(ii)+cc*tmp1(i,j)
          enddo
        enddo
      enddo
      deallocate(tmp1)
      end
