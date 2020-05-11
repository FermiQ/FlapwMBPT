      subroutine heg_vertex_p_psi
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
      if(psi_p=='00000') return
      if(psi_p=='60000') then
        call vertex_heg_p_only_spl(1,.false.)
c        call vertex_heg_spl(.true.,.false.,1,.false.)
      else if(psi_p=='66000') then
        call vertex_heg_p_only_spl(1,.true.)
      endif
      end