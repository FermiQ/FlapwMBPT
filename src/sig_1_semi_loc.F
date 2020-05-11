      subroutine sig_1_semi_loc(t_pw,grs,wsem,fif,s2,sig_semi,n,np)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: s2,n,np
      real*8, intent(in) :: grs(nrel,n,n,0:n_tau),fif(nrel,n,n,np),
     &                      wsem(np,np,ndim3_nu)
      complex*16, intent(in) :: t_pw(np,n,2,ndim3_tau)
      complex*16, intent(inout) :: sig_semi(n,n,2,ndim3_tau)
      integer :: it,i_tau,ind_nu,ind_tau,ii,j,i
      complex*16 :: cc
      complex*16, allocatable :: t_nu(:,:,:,:),wt(:,:,:,:),w(:,:),
     &                           wtt(:,:,:,:),gg(:,:),tt(:,:)
      allocate(t_nu(np,n,ndim3_nu,2))
      call tau_to_nu_full(t_nu,t_pw,np*n)
      allocate(wt(np,n,ndim3_nu,2))
      allocate(w(np,np))
      do it=1,2
        do ind_nu=1,ndim3_nu
          w=wsem(1:np,1:np,ind_nu)
          call zgemm('n','n',np,n,np,(1.d0,0.d0),w,np,
     &               t_nu(1,1,ind_nu,it),np,(0.d0,0.d0),
     &               wt(1,1,ind_nu,it),np)
        enddo
      enddo
      deallocate(w,t_nu)
      allocate(wtt(np,n,2,ndim3_tau))
      call nu_to_tau_cc(wt,wtt,np*n)
      deallocate(wt)
      allocate(gg(n,n))
      allocate(tt(np,n))
      do ind_tau=1,ndim3_tau
        do it=1,2
          i_tau=me_t*ndim3_tau+ind_tau-1
          if(it==2) i_tau=n_tau-i_tau
          if(irel/=2) then
            do j=1,n
              do i=1,n
                gg(i,j)=grs(1,i,j,i_tau)
              enddo
            enddo
          else if(irel==2) then
            do j=1,n
              do i=1,n
                gg(i,j)=dcmplx(grs(1,i,j,i_tau),grs(2,i,j,i_tau))
              enddo
            enddo
          endif
          call zgemm('n','t',np,n,n,(1.d0,0.d0),
     &               wtt(1,1,it,ind_tau),np,gg,n,(0.d0,0.d0),tt,np)
          do ii=1,n
            do j=1,n
              do i=1,np
                if(irel/=2) cc=fif(1,ii,j,i)
                if(irel==2) cc=dcmplx(fif(1,ii,j,i),fif(2,ii,j,i))
                sig_semi(ii,s2,it,ind_tau)=
     &            sig_semi(ii,s2,it,ind_tau)+cc*tt(i,j)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(wtt,gg,tt)
      end
