      subroutine k_sig0_stat(it,g_tau,s2,k1,ispin,k_sing)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,k1,s2,ispin
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2)
      complex*16, intent(out) :: k_sing(nbndf_bnd)
      integer :: ka,n,i
      real*8 :: w
      complex*16, allocatable :: gg(:,:)
      ka=k_a_from_c(k1)
      n=n_low_bnd(ka,ispin)
      w=8.d0*pi*q2aver_c/sqrt(amega)
c ------------------------------------------------------------------
      k_sing=(0.d0,0.d0)
      allocate(gg(n,n))
      call ferm_unpack_tau(gg,g_tau,n,n,nbndf_bnd,it)
      do i=1,n
        k_sing(i)=w*gg(i,s2)
      enddo
      deallocate(gg)
      end
