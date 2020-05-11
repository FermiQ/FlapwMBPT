      subroutine dphi_to_phi_large_r(e,b,s,alf,dalf,re,idm)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: idm
      real*8, intent(in) :: e,b,s,alf,dalf
      real*8, intent(out) :: re
      real*8 :: ada,lam,a
      if(idm==1) then
        re=-(s+b)*alf/s ! here alf - given LOG. derivative
      else if(idm==3) then
        if(irel/=2) then
          re=1.d0/s
          if(e<=0.d0) re=re+sqrt(-e)
          re=re*(s+b)
        else if(irel==2) then
          ada=dalf/2/alf
          a=ada**2-alf*e
          lam=-ada
          if(a>=0.d0) lam=lam+sqrt(a)
          re=lam*(s+b)+1.d0+b/s
        endif
      endif
      end