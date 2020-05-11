      subroutine metal_gw
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: nocc
      real*8 :: ho,lu
c ---------- Metal or Insulator ------------------------------------
      nocc=idnint(nelec)
      if(irel/=2.and.nspin==1.and.mod(nocc,2)==1) then
        metal=.true.
      else
        if(ubi==' hf') call hovb_lucb(e_bnd,ho,lu)
        if(ubi/=' hf') call hovb_lucb(e_qp,ho,lu)
        metal=.false.
        if(lu-ho<temperature*boltzman) metal=.true.
      endif
      if(maswrk) then
        if(metal) write(iun,*)
     &     ' The solid has been identified as METAL'
        if(.not.metal) write(iun,*)
     &     ' The solid has been identified as INSULATOR'
      endif
      end
