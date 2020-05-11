      subroutine v_coul_full_ii(nplwg,v_full,ind_q,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: nplwg,ind_q,khbr
      complex*16, intent(out) :: v_full(nplwg,nplwg)
      integer :: i
      v_full=(0.d0,0.d0)
      do i=1,nplwg
        if(khbr==0) v_full(i,i)=v_ii_3(i,ind_q)
        if(khbr==1) v_full(i,i)=v_ii_3_sr(i,ind_q)
      enddo
      end


      subroutine v_coul_full_ii_r(nplwg,v_full,ind_q,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: nplwg,ind_q,khbr
      real*8, intent(out) :: v_full(nplwg,nplwg)
      integer :: i
      v_full=0.d0
      do i=1,nplwg
        if(khbr==0) v_full(i,i)=v_ii_3(i,ind_q)
        if(khbr==1) v_full(i,i)=v_ii_3_sr(i,ind_q)
      enddo
      end
