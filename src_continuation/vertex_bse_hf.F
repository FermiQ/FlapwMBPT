      subroutine vertex_bse_hf
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      if(bse_kernel_p=='2000') then
c        call vertex_p_loc_hf(iter_ladder_p,0)
      else if(bse_kernel_p=='2200') then
c        call vertex_p_loc_hf(iter_ladder_p,iter_ladder_p)
      else if(bse_kernel_p=='4000') then
        call vertex_bnd_hf(iter_ladder_p,'wmt')
      else if(bse_kernel_p=='4400') then
        call vertex_bnd_hf(iter_ladder_p,'wmt')
      else if(bse_kernel_p=='6000') then
        call vertex_bnd_hf(iter_ladder_p,'bnd')
      else if(bse_kernel_p=='6600') then
        call vertex_bnd_hf(iter_ladder_p,'bnd')
      endif
      end
