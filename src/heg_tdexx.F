	  subroutine heg_alda(lda,chi_tmp)
	  use atom_mod
	  use etot_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      character*3, intent(in) :: chi_tmp
      integer, intent(in) :: lda
      real*8 :: th(4),pi2a,gro(3,3),gmod(3,3),romod(3),rolap(2)
      pi2a=(pi+pi)/par
      n_ixc_0=0
      if(nspin==1) then
        if(chi_tmp(1:1)=='1') then
          n_ixc_0=n_ixc_0+1
          ind_ixc(n_ixc_0)=1
        endif
        if(chi_tmp(3:3)=='1') then
          n_ixc_0=n_ixc_0+1
          ind_ixc(n_ixc_0)=2
        endif
      else if(nspin==2) then
        if(chi_tmp(1:1)=='1'.or.chi_tmp(3:3)=='1') then
          n_ixc_0=n_ixc_0+3
          ind_ixc(n_ixc_0-2)=1
          ind_ixc(n_ixc_0-1)=2
          ind_ixc(n_ixc_0)=3
        endif
        if(chi_tmp(2:2)=='1') then
          n_ixc_0=n_ixc_0+1
          ind_ixc(n_ixc_0)=4
        endif
      endif
c ------- Exchange-Correlation Kernel ------------------------------
      gro=0.d0
      gmod=0.d0
      romod=0.d0
      rolap=0.d0
      call theta_loc(lda,ro_spin(2),ro_spin(1),gro,gmod,romod,rolap,th,
     &               0)
      f_xc=f_xc+th(1)
      f_xc_zz=f_xc_zz+th(2)
      end
