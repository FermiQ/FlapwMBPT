      subroutine response_00_dft
      use atom_mod
      use etot_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use scalapack_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer :: ind_ch0(2),ind_ix0(3)
      call resp_lda_for_w
      if(nspin==1) then
        call resp_lda_w('___00',1,1,1,1)
      else if(nspin==2) then
        ind_ch0(1)=ind_chi(1)
        ind_ch0(2)=ind_chi(2) 
        ind_ix0(1)=1
        ind_ix0(2)=2
        ind_ix0(3)=3   
        call resp_lda_w('00_ZZ',2,ind_ch0,3,ind_ix0)
        if(lsclp) then
          deallocate(pz_sclp)
        else
          deallocate(p_0z_mm,p_0z_mi,p_0z_ii)
        endif
      endif
      deallocate(thet_me,theta_int)
      end
