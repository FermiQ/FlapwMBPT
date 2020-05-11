      real*8 function fermi_dirac(de)
	  use manager_mod
      implicit none
	  real*8 :: de,x,ffd
	  x=de*betta_t
	  if(x.lt.-100.d0) then
	    ffd=1.d0
	  else if(x.gt.100.d0) then
	    ffd=0.d0
	  else
	    ffd=1.d0/(1.d0+exp(x))
	  endif
	  fermi_dirac=ffd
      end


      real*8 function fermi_dirac_f(de)
      use manager_mod
      implicit none
      real*8 :: de,x,ffd
      x=de
      if(x.lt.-100.d0) then
        ffd=1.d0
      else if(x.gt.100.d0) then
        ffd=0.d0
      else
        ffd=1.d0/(1.d0+exp(x))
      endif
      fermi_dirac_f=ffd
      end


      real*8 function fermi_dirac_d(s,de)
      use manager_mod
      implicit none
      real*8 :: s,de,x,ffd
      x=de
      if(x.lt.-100.d0) then
        ffd=s
      else if(x.gt.100.d0) then
        ffd=0.d0
      else
        ffd=s/(1.d0+0.5d0*exp(x))
      endif
      fermi_dirac_d=ffd
      end
