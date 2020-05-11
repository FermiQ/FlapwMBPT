      subroutine get_r(ierror)
c	------- Mesh in the Unit Cell ----------------------------
	  use manager_mod
	  use parallel_mod
      use sitegen_mod
	  use solid_mod
      implicit none
	  integer, intent(out) :: ierror
      integer :: n1,n2,n3,ind,ind0,i1,i2,i3,ig,ir0,i,i0,n
      real*8 :: v(3),v1(3),v2(3),dx
      integer, allocatable :: num_s(:),ind_s(:,:),flag(:)
c  set divisions along x,y,z
      n1=nrdiv(1)-1
      n2=nrdiv(2)-1
      n3=nrdiv(3)-1
      if(allocated(ip_r)) deallocate(ip_r)
      allocate(ip_r(nr_full,ngroup))
      if(allocated(tshift_r)) deallocate(tshift_r)
      allocate(tshift_r(3,nr_full,ngroup))
      if(allocated(itshift_r)) deallocate(itshift_r)
      allocate(itshift_r(3,nr_full,ngroup))
      if(allocated(ir_ref)) deallocate(ir_ref)
      allocate(ir_ref(nr_full))
      if(allocated(r_group)) deallocate(r_group)
      allocate(r_group(nr_full))
      if(allocated(ind_r_star_back)) deallocate(ind_r_star_back)
      allocate(ind_r_star_back(nr_full))
c  initialize point-indicator array
	  ir_ref=0
c  *** find inequivalent r-points in the Unit Cell ***
      nr_red=0
c  do over whole IC
	  ind0=0
      allocate(num_s(nr_full))
      allocate(ind_s(nr_full,nr_full))
      num_s=0
      do i3=0,n3
        do i2=0,n2
          do i1=0,n1
	        ind0=ind0+1
	        r_pnt(1,ind0)=i1
	        r_pnt(2,ind0)=i2
	        r_pnt(3,ind0)=i3
            if(ir_ref(ind0)==0)then
              v=i1*r0b0(:,1)+i2*r0b0(:,2)+i3*r0b0(:,3)
c Reduce to zero cell ---------------------------------
c  rotate given point
              do ig=1,ngroup
                call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),2)
	            v1=v1+shift(:,ig)
c ----- Check the compatibility with symmetry operations -----
	            ierror=0
	            do ind=1,3
	              dx=dot_product(v1,q0b0(:,ind))
	              if(abs(nint(dx)-dx)>1.d-6) then
			        nrdiv(ind)=nrdiv(ind)+1
	                ierror=ierror+1
	              endif
	            enddo
	            if(ierror/=0) then
	              if(maswrk) write(iun,'(a35,3i5)')
     &	            ' Error in GET_R, NRDIV is set to be',
     &                            (nrdiv(ind),ind=1,3)
			      return
	            endif
	            call zone1_number(v1,q0b0,nrdiv,ind)
	            if(ir_ref(ind)==0) then
                  ir_ref(ind)=nr_red+1
	              r_group(ind)=ig
	            endif
              enddo
              nr_red=nr_red+1
              r_pnt_ind(nr_red)=ind0
              num_s(nr_red)=1
              ind_s(1,nr_red)=ind0
              ind_r_star_back(ind0)=1
            else
              ir0=ir_ref(ind0)
              num_s(ir0)=num_s(ir0)+1
              ind_s(num_s(ir0),ir0)=ind0
              ind_r_star_back(ind0)=num_s(ir0)
            endif
          enddo
        enddo
      enddo
c ------------- Calculation of IP_R and TSHIFT_R ----------
	  do ind0=1,nr_full
        v=r_pnt(1,ind0)*r0b0(:,1)+r_pnt(2,ind0)*r0b0(:,2)
     &			                 +r_pnt(3,ind0)*r0b0(:,3)
        do ig=1,ngroup
          call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),2)
	      v1=v1+shift(:,ig)
	      call zone1_number(v1,q0b0,nrdiv,ind)
          v2=r_pnt(1,ind)*r0b0(:,1)+r_pnt(2,ind)*r0b0(:,2)
     &			                   +r_pnt(3,ind)*r0b0(:,3)
	      ip_r(ind,ig)=ind0
	      tshift_r(:,ind,ig)=v2-v1
	      do i1=1,3
	        itshift_r(i1,ind,ig)=nint(dot_product(gbas(:,i1),
     &	                                          tshift_r(:,ind,ig)))
          enddo
        enddo
      enddo
c ------------- Information on the stars of r's ----------
      max_r_star=maxval(num_s)
      if(allocated(num_r_star)) deallocate(num_r_star)
      allocate(num_r_star(nr_red))
      if(allocated(ind_r_star)) deallocate(ind_r_star)
      allocate(ind_r_star(max_r_star,nr_red))
      do ir0=1,nr_red
        num_r_star(ir0)=num_s(ir0)
        do i1=1,num_s(ir0)
          ind_r_star(i1,ir0)=ind_s(i1,ir0)
        enddo
      enddo
      deallocate(num_s,ind_s)
      allocate(r_stars_based(nr_full))
      i=0
      do ir0=1,nr_red
        do i1=1,num_r_star(ir0)
          i=i+1
          r_stars_based(ind_r_star(i1,ir0))=i
        enddo
      enddo
      allocate(ip_r_star(max_r_star,nr_red,ngroup))
      do ir0=1,nr_red
        do i=1,num_r_star(ir0)
          ind=ind_r_star(i,ir0)
          do ig=1,ngroup
            i1=ip_r(ind,ig)
            ip_r_star(i,ir0,ig)=ind_r_star_back(i1)
          enddo
        enddo
      enddo
c ------------------------------------------------------------------
      if(.not.complex_ro) then
        allocate(num_r_inv(nr_red))
        allocate(ip_r_inv(max_r_star,nr_red))
        allocate(flag(max_r_star))
        nr_full_inv=0
        max_r_inv=0
        do ir0=1,nr_red
          n=0
          flag=0
          do i=1,num_r_star(ir0)
            if(flag(i)==0) then
              n=n+1
              flag(i)=1
              ip_r_inv(i,ir0)=n
              i0=ip_r_star(i,ir0,inv_num)
              if(i0/=i) then
                flag(i0)=1
                ip_r_inv(i0,ir0)=-n
              endif
            endif
          enddo
          nr_full_inv=nr_full_inv+n
          max_r_inv=max(max_r_inv,n)
          num_r_inv(ir0)=n
        enddo
        deallocate(flag)
        allocate(i_r_inv(nr_full_inv))
        allocate(ind_r_inv(max_r_inv,nr_red))
        n=0
        do ir0=1,nr_red
          do i=1,num_r_star(ir0)
            ind=ip_r_inv(i,ir0)
            if(ind>0) then
              n=n+1
              i_r_inv(n)=ind_r_star(i,ir0)
              ind_r_inv(ind,ir0)=n
            endif
          enddo
        enddo
        allocate(nr_full_ind(nr_full))
        do ir0=1,nr_red
          do i=1,num_r_star(ir0)
            n=ind_r_star(i,ir0)
            ind=ip_r_inv(i,ir0)
            i3=iabs(ind)
            i3=ind_r_inv(i3,ir0)
            if(ind>0) nr_full_ind(n)=i3
            if(ind<0) nr_full_ind(n)=-i3
          enddo
        enddo
      endif
      end
