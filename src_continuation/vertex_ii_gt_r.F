      subroutine vertex_ii_gt_r(i_tau,ll,t_rs,grs_ii,vx)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll
      real*8, intent(in) :: grs_ii(nr_red_red,nr_full_red,0:n_tau,
     &                             nqdiv_c)
      complex*16, intent(in) :: t_rs(nr_full_red,nr_full_red,nqdiv_c,ll)
      complex*16, intent(inout) :: vx(nr_full_red,nr_full_red,nqdiv_c,
     &                                ll)
      integer :: i,j,ir,ir0,ig,i0,j0
      real*8 :: v(3),v1(3)
      do i=1,nr_full_red
        i0=ir_ref_red(i)
        ig=ir_group_red(i)
        do ir=1,nqdiv_c
          do j=1,nr_full_red
            v=rr_coarse(:,ir)+tshift_r_red(:,i,ig)-tshift_r_red(:,j,ig)
            call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
            call zone1_number(v1,gbas,ndiv_c,ir0)
            j0=ip_r_red(j,ig)
            if(ll==1) then
              vx(i,j,ir,1)=vx(i,j,ir,1)
     &                    +grs_ii(i0,j0,i_tau,ir0)*t_rs(j,i,ir,1)
            else
              vx(i,j,ir,1)=vx(i,j,ir,1)
     &           +eqr(ir)*grs_ii(i0,j0,i_tau,ir0)*t_rs(j,i,ir,2)
              vx(i,j,ir,2)=vx(i,j,ir,2)
     &           +eqr(ir)*grs_ii(i0,j0,i_tau,ir0)*t_rs(j,i,ir,1)
            endif
          enddo  !! over j
        enddo   !! over ir
      enddo   !! over i
      end
