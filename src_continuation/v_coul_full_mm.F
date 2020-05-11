      subroutine v_coul_full_mm(v_full,ind_q,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ind_q,khbr
      complex*16, intent(out) :: v_full(n_pbmt,n_pbmt)
      integer :: ind_j,j,nn,n1,n2,i,ij
      complex*16, allocatable :: v0(:)
      nn=n_pbmt*(n_pbmt+1)/2
      n1=n0_tb_pbmt(me_t_b+1)+1
      n2=n1+nd_tb_pbmt(me_t_b+1)-1
      allocate(v0(nn))
      v0=(0.d0,0.d0)
      if(khbr==0) call zcopy(n2-n1+1,v_mm_3(1,ind_q),1,v0(n1),1)
      if(khbr==1) call zcopy(n2-n1+1,v_mm_3_sr(1,ind_q),1,v0(n1),1)
      if(nproc_t*nproc_b/=1) call dgop(v0,2*nn,'  +',comm_t_b)
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          v_full(i,j)=v0(ij)
          if(j/=i) v_full(j,i)=conjg(v0(ij))
        enddo
      enddo
      deallocate(v0)
      end


      subroutine v_coul_full_mm_r(v_full,ind_q,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ind_q,khbr
      real*8, intent(out) :: v_full(n_pbmt,n_pbmt)
      integer :: ind_j,j,nn,n1,n2,i,ij
      real*8, allocatable :: v0(:)
      nn=n_pbmt*(n_pbmt+1)/2
      n1=n0_tb_pbmt(me_t_b+1)+1
      n2=n1+nd_tb_pbmt(me_t_b+1)-1
      allocate(v0(nn))
      v0=0.d0
      if(khbr==0) call dcopy(n2-n1+1,v_mm_3(1,ind_q),1,v0(n1),1)
      if(khbr==1) call dcopy(n2-n1+1,v_mm_3_sr(1,ind_q),1,v0(n1),1)
      if(nproc_t*nproc_b/=1) call dgop(v0,nn,'  +',comm_t_b)
      ij=0
      do j=1,n_pbmt
        do i=1,j
          ij=ij+1
          v_full(i,j)=v0(ij)
          if(j/=i) v_full(j,i)=v0(ij)
        enddo
      enddo
      deallocate(v0)
      end
