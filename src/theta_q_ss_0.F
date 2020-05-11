      subroutine theta_q_ss_0(ka,n_pbt,thet_q)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n_pbt,ka
      real*8, intent(out) :: thet_q(ncmpl*n_pbt,n_pbt,n_ixc_0)
      integer :: iatom,isort,ind,n,j,i,j0,i0,gbs_number,ix,i1,j1
      real*8 :: v(3),v1(3)
      complex*16, allocatable :: tmp(:,:)
      thet_q=0.d0
c --------- MT part ------------------------------------------------
      do iatom=1,natom
        isort=is(iatom)
        ind=iopb(iatom)
        n=n_pbmt0(isort)
        do i=1,n_ixc_0
          if(ncmpl==1) then
            call dlacpy('a',n,n,thet_me(1,1,iatom,i),n_pbmtm,
     &                  thet_q(ind,ind,i),n_pbt)
          else
            allocate(tmp(n,n))
            tmp=thet_me(1:n,1:n,iatom,i)
            call zlacpy('a',n,n,tmp,n,thet_q(2*ind-1,ind,i),n_pbt)
            deallocate(tmp)
          endif
        enddo
      enddo   !!! over iatom
c --------- Interstitial -------------------------------------------
      n=nplwgw(ka)
      do j=1,n
        j0=indpw_gw(j,ka)  !! G' in GBS-list
        j1=n_pbmt+j
        do i=1,n
          i0=indpw_gw(i,ka)
          i1=n_pbmt+i
          v=gbs(:,j0)-gbs(:,i0)  !! G_i-G_j
          ind=gbs_number(v1)   !! G_i-G_j in GBS-list
          if(ncmpl==1) then
            do ix=1,n_ixc_0
              thet_q(i1,j1,ix)=theta_int(ind,ix)
            enddo
          else
            do ix=1,n_ixc_0
              call zcopy(1,theta_int(ind,ix),1,thet_q(2*i1-1,j1,ix),1)
            enddo
          endif
        enddo
      enddo
      end
