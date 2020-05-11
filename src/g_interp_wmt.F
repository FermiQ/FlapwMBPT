      subroutine g_interp_wmt(gexa,gspl_geom,gspl_asy,gasy,g_om_om,
     &                        g_nu_om,i1_omega,nom,nomnu,ispin,n,n1,
     &                        ir,indi,indj)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i1_omega,nom,nomnu,ispin,n,n1,indi,indj,ir
      real*8, intent(in) :: gexa(nfun_red,nfun_red,2,0:n_omega_exa,
     &                           nqdiv_c),
     &                      gspl_geom(nfun_red,nfun_red,2,
     &                                0:n_omega_geom+2,4,nqdiv_c),
     &                      gspl_asy(nfun_red,nfun_red,2,
     &                               n_omega_asy+1,4,nqdiv_c),
     &                      gasy(nfun_red,nfun_red,2,nqdiv_c)
      complex*16, intent(out) :: g_om_om(n,n1,2,nom),
     &                           g_nu_om(n,n1,2,nomnu)
      integer :: i_omega,i_nu,ir1,i,j
      real*8 :: w1,pib,v(3)
      complex*16, allocatable :: tmp(:,:)
      pib=pi/betta_t
      v=-rr_coarse(:,ir)
      call zone1_number(v,gbas,ndiv_c,ir1)
      allocate(tmp(n1,n))
      do i_omega=1,nom
        w1=w_om_adapt_omega(i_omega,i1_omega)
        call interp_g_wmt(gexa(1,1,1,0,ir),
     &                    gspl_geom(1,1,1,0,1,ir),
     &                    gspl_asy(1,1,1,1,1,ir),
     &                    gasy(1,1,1,ir),g_om_om(1,1,2,i_omega),
     &                    w1,n,n1,indi,indj,1)
        w1=-w_om_adapt_omega(i_omega,i1_omega)+w_omega(i1_omega)-pib
        if(w1>0.d0) then
          call interp_g_wmt(gexa(1,1,1,0,ir),
     &                      gspl_geom(1,1,1,0,1,ir),
     &                      gspl_asy(1,1,1,1,1,ir),
     &                      gasy(1,1,1,ir),
     &                      g_om_om(1,1,1,i_omega),w1,n,n1,indi,indj,1)
        else
          call interp_g_wmt(gexa(1,1,1,0,ir1),
     &                      gspl_geom(1,1,1,0,1,ir1),
     &                      gspl_asy(1,1,1,1,1,ir1),
     &                      gasy(1,1,1,ir1),
     &                      tmp,-w1,n1,n,indi,indj,1)
          do j=1,n1
            do i=1,n
              g_om_om(i,j,1,i_omega)=conjg(tmp(j,i))
            enddo
          enddo
        endif
      enddo  !! over i_omega
      do i_nu=1,nomnu
        w1=w_omega(i1_omega)-w_nu_adapt_omega(i_nu,i1_omega)
        if(w1>0.d0) then
          call interp_g_wmt(gexa(1,1,1,0,ir),
     &                      gspl_geom(1,1,1,0,1,ir),
     &                      gspl_asy(1,1,1,1,1,ir),
     &                      gasy(1,1,1,ir),
     &                      g_nu_om(1,1,2,i_nu),w1,n,n1,indi,indj,1)
        else
          call interp_g_wmt(gexa(1,1,1,0,ir1),
     &                      gspl_geom(1,1,1,0,1,ir1),
     &                      gspl_asy(1,1,1,1,1,ir1),
     &                      gasy(1,1,1,ir1),
     &                      tmp,-w1,n1,n,indi,indj,1)
          do j=1,n1
            do i=1,n
              g_nu_om(i,j,2,i_nu)=conjg(tmp(j,i))
            enddo
          enddo
        endif
        w1=w_nu_adapt_omega(i_nu,i1_omega)+pib
        call interp_g_wmt(gexa(1,1,1,0,ir),
     &                    gspl_geom(1,1,1,0,1,ir),
     &                    gspl_asy(1,1,1,1,1,ir),
     &                    gasy(1,1,1,ir),g_nu_om(1,1,1,i_nu),
     &                    w1,n,n1,indi,indj,1)
      enddo  !! over i_nu
      deallocate(tmp)
      end
