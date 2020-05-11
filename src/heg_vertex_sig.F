      subroutine heg_vertex_sig(mode,it_v)
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
	character*5, intent(in) :: mode
	integer, intent(in) :: it_v
c -------------- Vertex evaluation ----------------------------------
      sig_k_tau_vertex_heg=0.d0
      if(mode=='00000') return
      if(mode=='60000') then
        call vertex_heg_sig_1_1
c        call vertex_heg_spl(.false.,.true.,1,.false.)
      else if(mode=='66000') then
        call vertex_heg_spl(.false.,.true.,1,.true.)
      endif
      end