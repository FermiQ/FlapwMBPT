      subroutine heg_vertex_p_bse
	use atom_mod
	use etot_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
	implicit none
      p_q_nu_vertex_heg=0.d0
      if(bse_kernel_p=='0000') return
      if(bse_kernel_p=='6000') then
        call vertex_heg_p_only_spl(iter_ladder_p,.false.)
      else if(bse_kernel_p=='6600') then
        call vertex_heg_p_only_spl(iter_ladder_p,.true.)
      endif
      end