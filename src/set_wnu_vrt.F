      subroutine set_wnu_vrt
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: iq,ka,nb,ind_nu,i_nu,i,j
      complex*16, allocatable :: t(:,:)
      allocate(wnu_vrt(n_pbtot_red,n_pbtot_red,0:n_nu,npnt_c))
      wnu_vrt=0.d0
      do iq=1,npnt_c
        ka=k_a_from_c(iq)
        nb=n_pbmt_red+nplwgw_red(ka)
        allocate(t(nb,nb))
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          t=w_red_q(1:nb,1:nb,ind_nu,ka)
          if(iq==1) then
            do j=1,nb
              do i=1,nb
                t(i,j)=t(i,j)+pw_pb_tild_red(i,1)*ws_head_nu(i_nu)
     &                       *q2aver_c*conjg(pw_pb_tild_red(j,1))
              enddo
            enddo
          endif
          call pack_hermit(t,wnu_vrt(1,1,i_nu,iq),nb,nb,n_pbtot_red,
     &                     0.d0,1.d0)
        enddo  !! over ind_nu
        deallocate(t)
      enddo  !! over iq
      if(nproc_t/=1) call dgop(wnu_vrt,n_pbtot_red**2*(n_nu+1)*npnt_c,
     &                           '  +',comm_t)
      end
