      subroutine v_coul_full_mi(nplwg,v_full,ind_q,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: nplwg,ind_q,khbr
      complex*16, intent(out) :: v_full(n_pbmt,nplwg)
      integer :: n1,n2,j
      n1=nmpi_tb_n_pbmt(me_t_b+1)+1
      n2=n1+ndim_tb_n_pbmt(me_t_b+1)-1
      v_full=(0.d0,0.d0)
      if(ndim_tb_n_pbmt(me_t_b+1)/=0) then
        do j=1,nplwg
          if(khbr==0) call zcopy(n2-n1+1,v_mi_3(1,j,ind_q),1,
     &                           v_full(n1,j),1)
          if(khbr==1) call zcopy(n2-n1+1,v_mi_3_sr(1,j,ind_q),1,
     &                           v_full(n1,j),1)
        enddo
      endif
      if(nproc_t*nproc_b/=1) call dgop(v_full,2*n_pbmt*nplwg,'  +',
     &                                 comm_t_b)
      end


      subroutine v_coul_full_mi_r(nplwg,v_full,ind_q,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: nplwg,ind_q,khbr
      real*8, intent(out) :: v_full(n_pbmt,nplwg)
      integer :: n1,n2,j
      n1=nmpi_tb_n_pbmt(me_t_b+1)+1
      n2=n1+ndim_tb_n_pbmt(me_t_b+1)-1
      v_full=0.d0
      if(ndim_tb_n_pbmt(me_t_b+1)/=0) then
        do j=1,nplwg
          if(khbr==0) call dcopy(n2-n1+1,v_mi_3(1,j,ind_q),1,
     &                           v_full(n1,j),1)
          if(khbr==1) call dcopy(n2-n1+1,v_mi_3_sr(1,j,ind_q),1,
     &                           v_full(n1,j),1)
        enddo
      endif
      if(nproc_t*nproc_b/=1) call dgop(v_full,n_pbmt*nplwg,'  +',
     &                                 comm_t_b)
      end
