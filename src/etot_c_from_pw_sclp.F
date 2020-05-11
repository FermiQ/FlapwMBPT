      subroutine etot_c_from_pw_sclp(nrow,ncol,iq,i_nu,p,ws)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nrow,ncol,iq,i_nu
      complex*16, intent(in) :: p(nrow,ncol),ws(nrow,ncol)
      integer :: j,i
      real*8 :: de
      de=0.d0
      do j=1,ncol
        do i=1,nrow
          de=de-p(i,j)*conjg(ws(i,j))
        enddo
      enddo
      call dgop(de,1,'  +',comm_b)
      e_c_solid_pw=e_c_solid_pw+et_sum_nu(i_nu)*wgt(iq)*de
      end


      subroutine etot_c_from_pw_sclp_r(nrow,ncol,iq,i_nu,p,ws)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nrow,ncol,iq,i_nu
      real*8, intent(in) :: p(nrow,ncol),ws(nrow,ncol)
      integer :: j,i
      real*8 :: de
      de=0.d0
      do j=1,ncol
        do i=1,nrow
          de=de-p(i,j)*ws(i,j)
        enddo
      enddo
      call dgop(de,1,'  +',comm_b)
      e_c_solid_pw=e_c_solid_pw+et_sum_nu(i_nu)*wgt(iq)*de
      end
