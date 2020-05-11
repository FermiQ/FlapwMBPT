      subroutine heg_b_split
	use atom_mod
	use etot_mod
	use heg_mod
      use manager_mod
      use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      logical :: lb
      integer :: it
      b_split_heg=0.d0
      if(nspin==1) return
      if(b_extval>1.d-8) then
        b_split_heg=b_extval
        return
      endif
      lb=.false.
      if(ubi=='dft') then
        it=iter_h_ext/10000
        if(iter<=it) lb=.true.
      else if(ubi==' hf') then
        it=mod(iter_h_ext,10000)/1000
        if(iter-iter_dft<=it) lb=.true.
      else if(ubi==' gw') then
        it=mod(iter_h_ext,1000)/100
        if(iter-iter_dft-iter_hf<=it) lb=.true.
      else if(ubi==' qp') then
        it=mod(iter_h_ext,100)/10
        if(iter-iter_dft-iter_hf-iter_gw<=it) lb=.true.
      else if(ubi=='psi') then
        it=mod(iter_h_ext,10)
        if(iter-iter_dft-iter_hf-iter_gw-iter_qp<=it) lb=.true.
      else if(ubi=='bsp') then
        it=mod(iter_h_ext,10)
        if(iter-iter_dft-iter_hf-iter_gw-iter_qp-iter_psi<=it) lb=.true.
      endif
	if(.not.lb) return
	b_split_heg=magn_shift(1)
      end