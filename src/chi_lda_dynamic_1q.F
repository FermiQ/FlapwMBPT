      subroutine chi_lda_dynamic_1q(iq,n,nn,n_pbt,ch,n_ix0,ind_ix0,
     &                              thet_q,chi0,chi,wrx,nd,ns)
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
      integer, intent(in) :: n_ix0,ind_ix0(n_ix0),n,nn,n_pbt,iq,nd,ns
      real*8, intent(in) :: wrx(0:nrax_chi)
      complex*16, intent(in) :: thet_q(n_pbt,n_pbt,n_ixc_0)
      complex*16, intent(out) :: chi0(nn,nn,nd),chi(nn,nn,nd)
      integer :: i,info,j,ind_i
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
          if(chi_basis=='PB') then
            call chi0_lda_dynamic_1q(iq,2,1,chi0,n_pbt,wrx,nd,ns)
          else if(chi_basis=='PW') then
            call chi0_lda_dynamic_1qpw(iq,2,1,chi0,n_pbt,wrx,nd,ns)
          endif
          chi0=4.d0*chi0
        else if(ch=='___00'.or.ch=='___ZZ') then
          if(chi_basis=='PB') then
            call chi0_lda_dynamic_1q(iq,1,1,chi0,n_pbt,wrx,nd,ns)
          else if(chi_basis=='PW') then
            call chi0_lda_dynamic_1qpw(iq,1,1,chi0,n_pbt,wrx,nd,ns)
          endif
          chi0=2.d0*chi0
        endif
      else if(nn>n_pbt) then
      endif
c  ------- Enhanced Susceptibility ------------------------------
      do ind_i=1,nd
        i=ns+ind_i-1
        call zgemm('n','n',nn,nn,nn,(-1.d0,0.d0),chi0(1,1,ind_i),nn,
     &             ixc0,nn,(0.d0,0.d0),work,nn)
        do i=1,nn
          work(i,i)=(1.d0,0.d0)+work(i,i)
        enddo
        chi(:,:,ind_i)=chi0(:,:,ind_i)
        call zgesv(nn,nn,work,nn,ipiv,chi(1,1,ind_i),nn,info)
      enddo
      deallocate(ipiv,work,ixc0)
      end
