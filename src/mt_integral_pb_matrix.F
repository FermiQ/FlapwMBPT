      subroutine mt_integral_pb_matrix(mtp,pb,iq,chi_mt)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq
      real*8, intent(in) :: mtp(n_pbmt)
      complex*16, intent(in) :: pb(n_pbmt,n_pbmt)
      complex*16, intent(inout) :: chi_mt(natom)
      integer :: iatom,isort,n,ind,j,i
      complex*16 :: cc
      do iatom=1,natom
        isort=is(iatom)
        n=n_pbmt0(isort)-1
        ind=iopb(iatom)
        cc=(0.d0,0.d0)
        do j=ind,ind+n
          do i=ind,ind+n
            cc=cc+mtp(i)*pb(i,j)*mtp(j)
          enddo
        enddo
        chi_mt(iatom)=chi_mt(iatom)+wgt(iq)*cc
      enddo
      end
