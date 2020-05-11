      subroutine sig_1_semi_wmt(t_pw,grs,wnuq,fif,s2,ispin,sig_semi,
     &                          ir0,n,n1,np,indi,indj)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: s2,ispin,n,n1,np,ir0,indi,indj
      real*8, intent(in) :: grs(nrel,nfun_red,nfun_red,0:n_tau),
     &                      fif(nrel,maxel_red,maxel_red,n_pbmtm_red),
     &                      wnuq(n_pbmtm_red,n_pbmtm_red,ndim3_nu)
      complex*16, intent(in) :: t_pw(np,n1,2,ndim3_tau)
      complex*16, intent(inout) :: sig_semi(maxel_red,maxel_red,
     &                                      nrr_max_c,2,ndim3_tau)
      integer :: it,i_tau,ind_nu,ind_tau,ii,j,i
      complex*16 :: cc
      complex*16, allocatable :: t_nu(:,:,:,:),wt(:,:,:,:),w(:,:),
     &                           wtt(:,:,:,:),gg(:,:),tt(:,:)
      allocate(t_nu(np,n1,ndim3_nu,2))
      call tau_to_nu_full(t_nu,t_pw,np*n1)
      allocate(wt(np,n1,ndim3_nu,2))
      allocate(w(np,np))
      do it=1,2
        do ind_nu=1,ndim3_nu
          w=wnuq(1:np,1:np,ind_nu)
          call zgemm('n','n',np,n1,np,(1.d0,0.d0),w,np,
     &               t_nu(1,1,ind_nu,it),np,(0.d0,0.d0),
     &               wt(1,1,ind_nu,it),np)
        enddo
      enddo
      deallocate(w,t_nu)
      allocate(wtt(np,n1,2,ndim3_tau))
      call nu_to_tau_cc(wt,wtt,np*n1)
      deallocate(wt)
      allocate(gg(n,n1))
      allocate(tt(np,n))
      do ind_tau=1,ndim3_tau
        do it=1,2
          i_tau=me_t*ndim3_tau+ind_tau-1
          if(it==2) i_tau=n_tau-i_tau
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
     &                         grs(2,indi+i-1,indj+j-1,i_tau))
              enddo
            enddo
          endif
          call zgemm('n','t',np,n,n1,(1.d0,0.d0),
     &               wtt(1,1,it,ind_tau),np,gg,n,(0.d0,0.d0),tt,np)
          do ii=1,n
            do j=1,n
              do i=1,np
                if(irel/=2) cc=fif(1,ii,j,i)
                if(irel==2) cc=dcmplx(fif(1,ii,j,i),
     &                                fif(2,ii,j,i))
                sig_semi(ii,s2,ir0,it,ind_tau)=
     &            sig_semi(ii,s2,ir0,it,ind_tau)+cc*tt(i,j)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(wtt,gg,tt)
      end
