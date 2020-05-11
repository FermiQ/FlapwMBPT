      subroutine v_spectrum(iq,n0,ve,vie,ne)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq,n0
      integer, intent(out) :: ne
      real*8, intent(out) :: ve(n0)
      complex*16, intent(out) :: vie(n0,n0)
      integer :: i,j,npw,n,ind_q,iq0
      real*8 :: eps,pi2
      real*8, allocatable :: ve0(:)
      complex*16, allocatable :: vie0(:,:)
      npw=nplwgw(iq)
      n=n_pbmt+npw
      eps=eps_coul
      pi2=pi+pi
c ------ Full V matrix without divergent terms -------------------
      allocate(vie0(n,n))
      vie0=(0.d0,0.d0)
	  do ind_q=1,ndim3_k(me_k+1)
	    iq0=n3_mpi_k(me_k+1)+ind_q
        call v_coul_full(nplwgw(iq0),vie0,ind_q,0)
      enddo
	  if(nproc_k/=1) call dgop(vie0,2*n*n,'  +',comm_k)
      allocate(ve0(n))
      if(me_t_b==0) call eig_val_solver(n,n,vie0,ve0)
      if(nproc_t*nproc_b/=1) then
        call brdcst(msgdbl,vie0,16*n*n,0,comm_t_b)
        call brdcst(msgdbl,ve0,8*n,0,comm_t_b)
      endif
      j=0
      do i=1,n
          if(ve0(i)/ve0(n)>eps) cycle
          j=i
      enddo
      ne=n-j
      do i=1,ne
          ve(i)=ve0(n-i+1)
          vie(1:n,i)=vie0(:,n-i+1)
      enddo
      deallocate(vie0,ve0)
      end
