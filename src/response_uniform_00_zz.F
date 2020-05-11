      subroutine response_uniform_00_zz
	  use atom_mod
	  use etot_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      real*8 :: p0,p00,pzz,dpw_00,dpw_zz,xi
      call resp_uniform_prep_00_zz(p0)
      if(maswrk) then
        if(chi_cmp(1:1)/='0') write(iun,*)
     &   ' Non-interacting Uniform 00-polarizability is ',
     &              p0,' 1/Ry'
        if(chi_cmp(3:3)/='0') then
          write(iun,*)
     &     ' Non-interacting Uniform spin susceptibility is: '
          write(iun,*)' Volume: ',-3.34588d-4*p0
          xi=-2.98575d-5*p0*amega
          write(iun,*)'   Mass: ',xi/molar_mass,' cm^3/g'
          write(iun,*)'  Molar: ',xi,' cm^3/mol'
        endif
      endif
c --------- Leading term only vertex correction ---------------
      if(chi_cmp(1:1)/='0') call vertex_pol_uniform(0,dpw_00)
      if(chi_cmp(3:3)/='0') call vertex_pol_uniform(1,dpw_zz)
      if(chi_cmp(1:1)/='0') p00=p0+dpw_00
      if(chi_cmp(3:3)/='0') pzz=p0+dpw_zz
      if(maswrk) then
        if(chi_cmp(1:1)/='0') write(iun,*)
     &    ' Uniform 00-polarizability is ',p00,' 1/Ry'
        if(chi_cmp(3:3)/='0') then
          write(iun,*)
     &     ' Uniform spin susceptibility is: '
          write(iun,*)' Volume: ',-3.34588d-4*pzz
          xi=-2.98575d-5*pzz*amega
          write(iun,*)'   Mass: ',xi/molar_mass,' cm^3/g'
          write(iun,*)'  Molar: ',xi,' cm^3/mol'
        endif
      endif
      call timel('** RESPONSE_UNIFORM_00_ZZ finished *')
      end
