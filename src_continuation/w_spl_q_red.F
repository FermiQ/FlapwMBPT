      subroutine w_spl_q_red(iq,wnu,wexa,wspl_geom,wspl_asy,wasy)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq
      real*8, intent(in) :: wnu(n_pbtot_red,n_pbtot_red,0:n_nu,npnt_c)
	  real*8, intent(out) :: wexa(n_pbtot_red,n_pbtot_red,0:n_nu_exa),
     &                       wspl_geom(n_pbtot_red,n_pbtot_red,
     &                                 0:n_nu_geom+2,4),
     &                       wspl_asy(n_pbtot_red,n_pbtot_red,
     &                                n_nu_asy+1,4),
     &                       wasy(n_pbtot_red,n_pbtot_red)
	  integer :: ig,n,i_nu,ka,j,i,iq0
	  real*8, allocatable :: wnu0(:,:,:)
	  complex*16, allocatable :: tmp(:,:)
      wexa=0.d0
      wspl_geom=0.d0
      wspl_asy=0.d0
      wasy=0.d0
      iq0=i_kref_c(iq)
      ka=k_a_from_c(iq0)
      n=n_pbmt_red+nplwgw_red(ka)
	  ig=k_group_c(iq)
      allocate(wnu0(n,n,0:n_nu))
      allocate(tmp(n_pbtot_red,n_pbtot_red))
	  wnu0=0.d0
      do i_nu=0,n_nu
        call unpack_hermit(tmp,wnu(1,1,i_nu,iq0),n,n_pbtot_red,
     &                     n_pbtot_red,(0.d0,0.d0),(1.d0,0.d0))
	    call sym_w_red(iq,iq0,tmp)
        call pack_hermit(tmp,wnu0(1,1,i_nu),n,n_pbtot_red,n,0.d0,1.d0)
      enddo  !! over i_nu
      deallocate(tmp)
      do j=1,n
        do i=1,n
          call spline_inhmg(w_nu(n_nu_exa-1),wnu0(i,j,n_nu_exa-1:),
     &                      wspl_geom(i,j,:,1),wspl_geom(i,j,:,2),
     &                      wspl_geom(i,j,:,3),wspl_geom(i,j,:,4),
     &                      n_nu_geom+3,0,0.d0,0.d0)
          call spline_inhmg(xm_nu(0),wnu0(i,j,n_nu_exa+n_nu_geom:),
     &                      wspl_asy(i,j,:,1),wspl_asy(i,j,:,2),
     &                      wspl_asy(i,j,:,3),wspl_asy(i,j,:,4),
     &                      n_nu_asy+1,0,0.d0,0.d0)
          wexa(i,j,:)=wnu0(i,j,0:n_nu_exa)
          wasy(i,j)=wnu0(i,j,n_nu)
        enddo
      enddo 
      deallocate(wnu0)
      end
