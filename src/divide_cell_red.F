      subroutine divide_cell_red
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer :: m,ind,i3,i2,i1,ind0,ig
	real*8 :: rad_inscribed_sphere,v(3),v1(3),v2(3)
	nr_full_red=nrdiv_red(1)*nrdiv_red(2)*nrdiv_red(3)
      r0b0_red(:,1) = rbas(:,1)/nrdiv_red(1)
      r0b0_red(:,2) = rbas(:,2)/nrdiv_red(2)
      r0b0_red(:,3) = rbas(:,3)/nrdiv_red(3)
      q0b0_red(:,1) = gbas(:,1)*nrdiv_red(1)
      q0b0_red(:,2) = gbas(:,2)*nrdiv_red(2)
      q0b0_red(:,3) = gbas(:,3)*nrdiv_red(3)
      allocate(r_pnt_red(3,nr_full_red))
      ind=0
      do i3=0,nrdiv_red(3)-1
        do i2=0,nrdiv_red(2)-1
          do i1=0,nrdiv_red(1)-1
	      ind=ind+1
	      r_pnt_red(1,ind)=i1
	      r_pnt_red(2,ind)=i2
	      r_pnt_red(3,ind)=i3
	    enddo
	  enddo
	enddo
	cut_pw_red=rad_inscribed_sphere(q0b0_red)
      cut_bnd_red=cut_pw_red*cut_vrt_ratio
c ---------------------------------------------------------------------
	if(maswrk) then
	  write(iun,'(a11,3i4)')'NRDIV_RED: ',(nrdiv_red(m),m=1,3)
	  write(iun,'(a13,i5)')'NR_FULL_RED: ',nr_full_red
	  write(iun,'(a13,f12.5)')'CUT_BND_RED: ',cut_bnd_red
	  write(iun,'(a12,f12.5)')'CUT_PW_RED: ',cut_pw_red
	endif
c ------------- Calculation of IP_R_RED and TSHIFT_R_RED ----------
      if(allocated(ip_r_red)) deallocate(ip_r_red)
      allocate(ip_r_red(nr_full_red,ngroup))
      if(allocated(tshift_r_red)) deallocate(tshift_r_red)
      allocate(tshift_r_red(3,nr_full_red,ngroup))
	do ind0=1,nr_full_red
        v=r_pnt_red(1,ind0)*r0b0_red(:,1)
     &   +r_pnt_red(2,ind0)*r0b0_red(:,2)
     &   +r_pnt_red(3,ind0)*r0b0_red(:,3)
        do ig=1,ngroup
          call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),2)
	    v1=v1+shift(:,ig)
	    call zone1_number(v1,q0b0_red,nrdiv_red,ind)
          v2=r_pnt_red(1,ind)*r0b0_red(:,1)
     &      +r_pnt_red(2,ind)*r0b0_red(:,2)
     &      +r_pnt_red(3,ind)*r0b0_red(:,3)
	    ip_r_red(ind,ig)=ind0
	    tshift_r_red(:,ind,ig)=v2-v1
        enddo
      enddo
c  initialize point-indicator array
      allocate(ir_ref_red(nr_full_red))
	ir_ref_red=0
      allocate(ir_group_red(nr_full_red))
	ir_group_red=0
c  *** find inequivalent r-points in the Unit Cell ***
      nr_red_red=0
      allocate(r_pnt_ind_red(nr_full_red))
	r_pnt_ind_red=0
	do ind0=1,nr_full_red
        if(ir_ref_red(ind0)==0)then
          v=r_pnt_red(1,ind0)*r0b0_red(:,1)
     &     +r_pnt_red(2,ind0)*r0b0_red(:,2)
     &     +r_pnt_red(3,ind0)*r0b0_red(:,3)
c  rotate given point
          do ig=1,ngroup
            call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),2)
	      v1=v1+shift(:,ig)
	      call zone1_number(v1,q0b0_red,nrdiv_red,ind)
	      if(ir_ref_red(ind)==0) then
              ir_ref_red(ind)=nr_red_red+1
	        ir_group_red(ind)=ig
	      endif
          enddo
          nr_red_red=nr_red_red+1
          r_pnt_ind_red(nr_red_red)=ind0
        endif
      enddo
      call timel('*** Divide_CELL_RED   finished *****')
      end