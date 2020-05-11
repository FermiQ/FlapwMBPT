      subroutine w_spl_loc_red(iatom,wnu,wexa,wspl_geom,wspl_asy,wasy)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iatom
      real*8, intent(in) :: wnu(n_pbmtm_red,n_pbmtm_red,0:n_nu)
	  real*8, intent(out) :: wexa(n_pbmtm_red,n_pbmtm_red,0:n_nu_exa),
     &                       wspl_geom(n_pbmtm_red,n_pbmtm_red,
     &                                 0:n_nu_geom+2,4),
     &                       wspl_asy(n_pbmtm_red,n_pbmtm_red,
     &                                n_nu_asy+1,4),
     &                       wasy(n_pbmtm_red,n_pbmtm_red)
	  integer :: n,j,i,isort
      wexa=0.d0
      wspl_geom=0.d0
      wspl_asy=0.d0
      wasy=0.d0
      isort=is(iatom)
      n=n_pbmt0_red(isort)
      do j=1,n
        do i=1,n
          call spline_inhmg(w_nu(n_nu_exa-1),wnu(i,j,n_nu_exa-1:),
     &                      wspl_geom(i,j,:,1),wspl_geom(i,j,:,2),
     &                      wspl_geom(i,j,:,3),wspl_geom(i,j,:,4),
     &                      n_nu_geom+3,0,0.d0,0.d0)
          call spline_inhmg(xm_nu(0),wnu(i,j,n_nu_exa+n_nu_geom:),
     &                      wspl_asy(i,j,:,1),wspl_asy(i,j,:,2),
     &                      wspl_asy(i,j,:,3),wspl_asy(i,j,:,4),
     &                      n_nu_asy+1,0,0.d0,0.d0)
          wexa(i,j,:)=wnu(i,j,0:n_nu_exa)
          wasy(i,j)=wnu(i,j,n_nu)
        enddo
      enddo 
      end
