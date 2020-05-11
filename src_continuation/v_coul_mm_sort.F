      subroutine v_coul_mm_sort(v_full,ind_q,v_mm,n1,n2,ind0,jnd0)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ind_q,n1,n2,ind0,jnd0
      complex*16, intent(in) :: v_mm(nd_tb_pbmt(me_t_b+1))
      complex*16, intent(out) :: v_full(n1,n2)
      integer :: ij,j,i,n00,n11,jj
      complex*16, allocatable :: v0(:,:)
      n00=n0_tb_pbmt(me_t_b+1)
      n11=n00+nd_tb_pbmt(me_t_b+1)
      allocate(v0(n1,n2))
      v0=(0.d0,0.d0)
      if(ind0<=jnd0) then
        do j=jnd0+1,jnd0+n2
          jj=j*(j-1)/2
          do i=ind0+1,min(ind0+n1,j)
            ij=jj+i
            if(ij<=n00.or.ij>n11) cycle
            v0(i-ind0,j-jnd0)=v_mm(ij-n00)
          enddo
        enddo
        if(ind0==jnd0) then
          do j=1,n2
            do i=1,j-1
              v0(j,i)=conjg(v0(i,j))
            enddo
          enddo
        endif
      else
        do j=ind0+1,ind0+n1
          jj=j*(j-1)/2
          do i=jnd0+1,jnd0+n2
            ij=jj+i
            if(ij<=n00.or.ij>n11) cycle
            v0(j-ind0,i-jnd0)=conjg(v_mm(ij-n00))
          enddo
        enddo
      endif
      if(nproc_t*nproc_b/=1) call dgop(v0,2*n1*n2,'  +',comm_t_b)
      v_full=v0
      deallocate(v0)
      end
