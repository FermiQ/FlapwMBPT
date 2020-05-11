      subroutine g_interp_bnd(gexa,gspl_geom,gspl_asy,gasy,g_om_nu,i_nu,
     &                        nom)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i_nu,nom
	  real*8, intent(in) :: gexa(nbndf_bnd,nbndf_bnd,2,0:n_omega_exa,
     &                           npnt_c,nspin),
     &                      gspl_geom(nbndf_bnd,nbndf_bnd,2,
     &                                0:n_omega_geom+2,4,npnt_c,nspin),
     &                      gspl_asy(nbndf_bnd,nbndf_bnd,2,
     &                               n_omega_asy+1,4,npnt_c,nspin),
     &                      gasy(nbndf_bnd,nbndf_bnd,2,npnt_c,nspin)
      complex*16, intent(out) :: g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,
     &                                   npnt_c,nspin)
      integer :: k,ispin,nb,i_omega,kf
      real*8 :: w1
c ------ G-interpolation ------------------------------------------
      do ispin=1,nspin
        do k=1,npnt_c
          kf=k_a_from_c(k)
	    nb=n_low_bnd(kf,ispin)
	    do i_omega=1,nom
            w1=w_om_adapt_nu(i_omega,i_nu)
            call interp_g_sq_spl(gexa(1,1,1,0,k,ispin),
     &                           gspl_geom(1,1,1,0,1,k,ispin),
     &                           gspl_asy(1,1,1,1,1,k,ispin),
     &                           gasy(1,1,1,k,ispin),
     &                           g_om_nu(1,1,1,i_omega,k,ispin),w1,
     &                           nbndf_bnd,1)
            w1=w_om_adapt_nu(i_omega,i_nu)-w_nu(i_nu)
            call interp_g_sq_spl(gexa(1,1,1,0,k,ispin),
     &                           gspl_geom(1,1,1,0,1,k,ispin),
     &                           gspl_asy(1,1,1,1,1,k,ispin),
     &                           gasy(1,1,1,k,ispin),
     &                           g_om_nu(1,1,2,i_omega,k,ispin),w1,
     &                           nbndf_bnd,1)
          enddo  !! over ind_omega
        enddo  !! over k
      enddo   !! over ispin
      end
      
      
      
      
      subroutine g_interp_sig(gexa,gspl_geom,gspl_asy,gasy,g_om_om,
     &                        g_nu_om,i1_omega,nom,nomnu,ispin)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i1_omega,nom,nomnu,ispin
	  real*8, intent(in) :: gexa(nbndf_bnd,nbndf_bnd,2,0:n_omega_exa,
     &                           npnt_c),
     &                      gspl_geom(nbndf_bnd,nbndf_bnd,2,
     &                                0:n_omega_geom+2,4,npnt_c),
     &                      gspl_asy(nbndf_bnd,nbndf_bnd,2,
     &                               n_omega_asy+1,4,npnt_c),
     &                      gasy(nbndf_bnd,nbndf_bnd,2,npnt_c)
      complex*16, intent(out) :: g_om_om(nbndf_bnd,nbndf_bnd,2,nom,
     &                                   npnt_c),
     &                           g_nu_om(nbndf_bnd,nbndf_bnd,2,nomnu,
     &                                   npnt_c)
      integer :: k,nb,i_omega,kf,i_nu
      real*8 :: w1,pib
      pib=pi/betta_t
c ------ G-interpolation ------------------------------------------
      do k=1,npnt_c
        kf=k_a_from_c(k)
	    nb=n_low_bnd(kf,ispin)
	    do i_omega=1,nom
          w1=w_om_adapt_omega(i_omega,i1_omega)
          call interp_g_sq_spl(gexa(1,1,1,0,k),gspl_geom(1,1,1,0,1,k),
     &                         gspl_asy(1,1,1,1,1,k),gasy(1,1,1,k),
     &                         g_om_om(1,1,2,i_omega,k),w1,nbndf_bnd,1)
          w1=-w_om_adapt_omega(i_omega,i1_omega)+w_omega(i1_omega)-pib
          call interp_g_sq_spl(gexa(1,1,1,0,k),gspl_geom(1,1,1,0,1,k),
     &                         gspl_asy(1,1,1,1,1,k),gasy(1,1,1,k),
     &                         g_om_om(1,1,1,i_omega,k),w1,nbndf_bnd,1)
        enddo  !! over i1_omega
	    do i_nu=1,nomnu
          w1=w_omega(i1_omega)-w_nu_adapt_omega(i_nu,i1_omega)
          call interp_g_sq_spl(gexa(1,1,1,0,k),gspl_geom(1,1,1,0,1,k),
     &                         gspl_asy(1,1,1,1,1,k),gasy(1,1,1,k),
     &                         g_nu_om(1,1,2,i_nu,k),w1,nbndf_bnd,1)
          w1=w_nu_adapt_omega(i_nu,i1_omega)+pib
          call interp_g_sq_spl(gexa(1,1,1,0,k),gspl_geom(1,1,1,0,1,k),
     &                         gspl_asy(1,1,1,1,1,k),gasy(1,1,1,k),
     &                         g_nu_om(1,1,1,i_nu,k),w1,nbndf_bnd,1)
        enddo  !! over i_nu
      enddo  !! over k
      end
