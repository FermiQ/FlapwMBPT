      subroutine tos_trnsp(x,ttos,eb,n,beta,np,wg,e_don,s_don)
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,np
      real*8, intent(in) :: eb(n,np,nspin),beta,wg(np),e_don,s_don
      real*8, intent(out) :: ttos
      integer :: ispin,k,ibnd
      real*8 :: const,st,fermi_dirac_f,x,fermi_dirac_d
      const=2.d0/dfloat(nspin*nrel)
      ttos=0.d0
      do ispin=1,nspin
        do k=1,np
          st=0.d0
          do ibnd=1,n
            st=st+fermi_dirac_f((eb(ibnd,k,ispin)-x)*beta)
          enddo
          ttos=ttos+wg(k)*st*const
        enddo  !!! over k
      enddo  !!! over ispin
      ttos=ttos+fermi_dirac_d(s_don,(e_don-x)*beta)
      end
