      subroutine vertex_ii_gt_c(i_tau,ll,t_rs,grs_ii,vx)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll
      complex*16, intent(in) :: grs_ii(2*nr_red_red,2*nr_full_red,
     &                                 0:n_tau,nqdiv_c)
      complex*16, intent(in) :: t_rs(nr_full_red,nr_full_red,nqdiv_c,ll)
      complex*16, intent(inout) :: vx(2*nr_full_red,2*nr_full_red,
     &                                nqdiv_c,ll)
      integer :: i,j,j1,i1,ir,ir0,ig,i0,j0,irl,jrl,is1,js,is0,js0
      real*8 :: v(3),v1(3)
      complex*16 :: cc
      do i=1,nr_full_red
        i0=ir_ref_red(i)
        ig=ir_group_red(i)
        do ir=1,nqdiv_c
          do j=1,nr_full_red
            v=rr_coarse(:,ir)+tshift_r_red(:,i,ig)-tshift_r_red(:,j,ig)
            call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
            call zone1_number(v1,gbas,ndiv_c,ir0)
            j0=ip_r_red(j,ig)
            do jrl=1,2
              j1=(jrl-1)*nr_full_red+j
              do irl=1,2
                i1=(irl-1)*nr_full_red+i
                do js=1,2
                  js0=(js-1)*nr_full_red+j0
                  do is1=1,2
                    is0=(is1-1)*nr_red_red+i0
                    cc=conjg(uj((is1-1)*2+irl,ig))
     &                      *uj((js-1)*2+jrl,ig)
                    if(ll==1) then
                      vx(i1,j1,ir,1)=vx(i1,j1,ir,1)
     &                      +grs_ii(is0,js0,i_tau,ir0)*t_rs(j,i,ir,1)
                    else
                      vx(i1,j1,ir,1)=vx(i1,j1,ir,1)
     &                 +eqr(ir)*grs_ii(is0,js0,i_tau,ir0)*t_rs(j,i,ir,2)
                      vx(i1,j1,ir,2)=vx(i1,j1,ir,2)
     &                 +eqr(ir)*grs_ii(is0,js0,i_tau,ir0)*t_rs(j,i,ir,1)
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo  !! over j
        enddo   !! over ir
      enddo   !! over i
      end
