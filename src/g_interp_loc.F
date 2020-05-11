      subroutine g_interp_loc(gexa,gspl_geom,gspl_asy,gasy,g_om_nu,i_nu,
     &                        nom,n)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_nu,nom,n
      real*8, intent(in) :: gexa(n,n,2,0:n_omega_exa,nspin),
     &                      gspl_geom(n,n,2,0:n_omega_geom+2,4,nspin),
     &                      gspl_asy(n,n,2,n_omega_asy+1,4,nspin),
     &                      gasy(n,n,2,nspin)
      complex*16, intent(out) :: g_om_nu(n,n,2,nom,nspin)
      integer :: ispin,i_omega,j,i
      real*8 :: w1
      complex*16, allocatable :: tmp(:,:)
      allocate(tmp(n,n))
c ------ G-interpolation ------------------------------------------
      do ispin=1,nspin
        do i_omega=1,nom
          w1=w_om_adapt_nu(i_omega,i_nu)
          call interp_g_loc(gexa(1,1,1,0,ispin),
     &                      gspl_geom(1,1,1,0,1,ispin),
     &                      gspl_asy(1,1,1,1,1,ispin),
     &                      gasy(1,1,1,ispin),
     &                      g_om_nu(1,1,1,i_omega,ispin),w1,
     &                      n,1)
          w1=w_om_adapt_nu(i_omega,i_nu)-w_nu(i_nu)
          if(w1>0.d0) then
            call interp_g_loc(gexa(1,1,1,0,ispin),
     &                        gspl_geom(1,1,1,0,1,ispin),
     &                        gspl_asy(1,1,1,1,1,ispin),
     &                        gasy(1,1,1,ispin),
     &                        g_om_nu(1,1,2,i_omega,ispin),w1,
     &                        n,1)
          else
            call interp_g_loc(gexa(1,1,1,0,ispin),
     &                        gspl_geom(1,1,1,0,1,ispin),
     &                        gspl_asy(1,1,1,1,1,ispin),
     &                        gasy(1,1,1,ispin),tmp,-w1,n,1)
            do j=1,n
              do i=1,n
                g_om_nu(i,j,2,i_omega,ispin)=conjg(tmp(j,i))
              enddo
            enddo
        endif
        enddo  !! over ind_omega
      enddo   !! over ispin
      deallocate(tmp)
      end
