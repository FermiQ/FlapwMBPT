      subroutine chi_lda_static_1q(iq,n,nn,n_pbt,ch,n_ix0,ind_ix0,
     &                             thet_q,chi0,chi)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ch
      integer, intent(in) :: n_ix0,ind_ix0(n_ix0),n,nn,n_pbt,iq
      complex*16, intent(in) :: thet_q(n_pbt,n_pbt,n_ixc_0)
      complex*16, intent(out) :: chi0(nn,nn),chi(nn,nn)
      integer :: i,info,j
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: work(:,:),ixc0(:,:)
      allocate(ipiv(nn))
      allocate(work(nn,nn))
      allocate(ixc0(nn,nn))
c -------------------------------------------------------------------
      if(n==1.or.ch=='___XY') then
        j=ind_ix0(1)
        ixc0(:,:)=thet_q(:,:,j)
        if(ch=='___XY') ixc0(:,:)=ixc0(:,:)/2 
      else
        j=ind_ix0(1)
        ixc0(1:nn/2,1:nn/2)=thet_q(:,:,j)
        j=ind_ix0(2)
        ixc0(nn/2+1:nn,nn/2+1:nn)=thet_q(:,:,j)
        j=ind_ix0(3)
        ixc0(nn/2+1:nn,1:nn/2)=thet_q(:,:,j)
        ixc0(1:nn/2,nn/2+1:nn)=thet_q(:,:,j)
      endif
c --------------------------------------------------------------------
      if(nn==n_pbt) then
        if(ch=='___XY') then
          call chi0_lda_static_1q(iq,1,2,chi0,n_pbt)
          chi0=4.d0*chi0
        else if(ch=='___00'.or.ch=='___ZZ') then
          call chi0_lda_static_1q(iq,1,1,chi0,n_pbt)
          chi0=2.d0*chi0
        endif
      else if(nn>n_pbt) then
      endif
c  ------- Enhanced Susceptibility ------------------------------
      call zgemm('n','n',nn,nn,nn,(-1.d0,0.d0),chi0,nn,ixc0,nn,
     &           (0.d0,0.d0),work,nn)
      do i=1,nn
        work(i,i)=(1.d0,0.d0)+work(i,i)
      enddo
      chi=chi0
      call zgesv(nn,nn,work,nn,ipiv,chi,nn,info)
      deallocate(ipiv,work,ixc0)
      end