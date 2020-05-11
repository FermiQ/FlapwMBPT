      subroutine tsig0_end_stat(t_sing,gctau,sig,nk)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nk
      complex*16, intent(in) :: t_sing(nbndf_bnd),
     &                          gctau(nbndf_bnd,nbndf_bnd)
      complex*16, intent(inout) :: sig(nbndf_bnd)
      integer :: i,j
      real*8 :: w
      complex*16, allocatable :: tt(:)
      w=8.d0*pi*q2aver_c
      allocate(tt(nk))
      tt=(0.d0,0.d0)
      do j=1,nk
        do i=1,nk
          tt(i)=tt(i)+gctau(i,j)*t_sing(j)
        enddo
      enddo
      tt=tt*w
      w=1.d0/sqrt(amega)/nqdiv_c
      do i=1,nk
        sig(i)=sig(i)+w*tt(i)
      enddo
      deallocate(tt)
      end
