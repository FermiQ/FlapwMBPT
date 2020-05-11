      subroutine sig_1_stat_wmt(it,ind_tau,t_pw,vloc,fif,s2,
     &                          ispin,grs,sig_stat,n,n1,np,np1,iatom,
     &                          jatom,ir)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,ind_tau,s2,ispin,n,n1,np,np1,iatom,
     &                       jatom,ir
      real*8, intent(in) :: grs(nrel,nfun_red,nfun_red,0:n_tau),
     &                      fif(nrel,maxel_red,maxel_red,n_pbmtm_red,
     &                          nsort),
     &                      vloc(n_pbmtm_red,n_pbmtm_red,natom)
      complex*16, intent(inout) :: sig_stat(n)
      complex*16, intent(out) :: t_pw(np,n1)
      integer :: i_tau,i1_tau,i,ii,indi,indj,jsort,isort,j
      complex*16 :: cc
      complex*16, allocatable :: k_pw(:,:),w(:,:),tmp(:,:),tmp1(:,:),
     &                           gg(:,:)
      isort=is(iatom)
      jsort=is(jatom)
      indi=io_lem_red(iatom)
      indj=io_lem_red(jatom)
      i_tau=me_t*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau
      i1_tau=n_tau-i_tau
      allocate(k_pw(n,np1))
      call k_sig_stat_wmt(it,ind_tau,grs(1,1,1,i_tau),
     &                    fif(1,1,1,1,jsort),vloc(1,1,jatom),s2,ispin,
     &                    k_pw,n,n1,np1,iatom,jatom)
      if(irel/=2) then
        call tsig_wmt(grs,k_pw,t_pw,fif,i1_tau,n,n1,np,np1,iatom,
     &                jatom)
      else
        call tsig_wmt_c(grs,k_pw,t_pw,fif,i1_tau,n,n1,np,np1,iatom,
     &                  jatom)
      endif
      deallocate(k_pw)
      allocate(tmp(np,n1))
      allocate(tmp1(np,n))
      allocate(w(np,np))
      w=vloc(1:np,1:np,iatom)
      call zgemm('n','n',np,n1,np,(1.d0,0.d0),w,np,t_pw,np,(0.d0,0.d0),
     &           tmp,np)
      deallocate(w)
      allocate(gg(n,n1))
      if(irel/=2) then
        do j=1,n1
          do i=1,n
            gg(i,j)=grs(1,indi+i-1,indj+j-1,i_tau)
          enddo
        enddo
      else if(irel==2) then
        do j=1,n1
          do i=1,n
            gg(i,j)=dcmplx(grs(1,indi+i-1,indj+j-1,i_tau),
     &                     grs(2,indi+i-1,indj+j-1,i_tau))
          enddo
        enddo
      endif
      call zgemm('n','t',np,n,n1,(1.d0,0.d0),tmp,np,gg,n,(0.d0,0.d0),
     &           tmp1,np)
      deallocate(gg,tmp)
      do ii=1,n
        do j=1,n
          do i=1,np
            if(irel/=2) cc=fif(1,ii,j,i,isort)
            if(irel==2) cc=dcmplx(fif(1,ii,j,i,isort),
     &                            fif(2,ii,j,i,isort))
            sig_stat(ii)=sig_stat(ii)+cc*tmp1(i,j)
          enddo
        enddo
      enddo
      deallocate(tmp1)
      end
