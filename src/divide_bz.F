      subroutine divide_bz
	  use atom_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer :: m,ierror
	  real*8 :: rad_inscribed_sphere,qb(3,3),ct0
2	  nr_full=nrdiv(1)*nrdiv(2)*nrdiv(3)
      r0b0(:,1) = rbas(:,1)/nrdiv(1)
      r0b0(:,2) = rbas(:,2)/nrdiv(2)
      r0b0(:,3) = rbas(:,3)/nrdiv(3)
      q0b0(:,1) = gbas(:,1)*nrdiv(1)
      q0b0(:,2) = gbas(:,2)*nrdiv(2)
      q0b0(:,3) = gbas(:,3)*nrdiv(3)
	  if(allocated(r_pnt)) deallocate(r_pnt)
	  if(allocated(r_pnt_ind)) deallocate(r_pnt_ind)
      allocate(r_pnt_ind(nr_full))
      allocate(r_pnt(3,nr_full))
c ------------ Here we set up r-mesh for given NRDIV     ----------
	  call get_r(ierror)
	  if(ierror/=0) goto 2
c --------- Setting up CUTOFF and CUTOFFGW for given NRDIV  ---
      qb(:,1) = gbas(:,1)*mdiv(1)
      qb(:,2) = gbas(:,2)*mdiv(2)
      qb(:,3) = gbas(:,3)*mdiv(3)
	  cutoffro=rad_inscribed_sphere(qb)
	  ct0=rad_inscribed_sphere(q0b0)
	  cutoffgw=ct0
	  cutoff=cut_gw_ratio*ct0
	  if(maswrk) then
	    write(iun,'(a10,f8.4,a12,f8.4,a12,f8.4)')' CUTOFF = ',cutoff,
     &	               ' CUTOFFRO = ',cutoffro,' CUTOFFGW = ',cutoffgw
      endif
c --------------------------------------------------------------
12	  nqdiv=ndiv(1)*ndiv(2)*ndiv(3)
	  if(allocated(pnt)) deallocate(pnt)
	  if(allocated(wgt)) deallocate(wgt)
      allocate(pnt(3,nqdiv))
      allocate(wgt(nqdiv))
	  do m = 1, 3
        rb0(m,1) = rbas(m,1)*ndiv(1)
        rb0(m,2) = rbas(m,2)*ndiv(2)
        rb0(m,3) = rbas(m,3)*ndiv(3)
        qb0(m,1) = gbas(m,1)/ndiv(1)
        qb0(m,2) = gbas(m,2)/ndiv(2)
        qb0(m,3) = gbas(m,3)/ndiv(3)
      enddo
      call get_k(ierror)
	  if(ierror/=0) then
	    ndiv(ierror)=ndiv(ierror)+1
	    if(maswrk) write(iun,'(a34,3i5)')
     &	            ' Error in GET_K, NDIV is set to be',
     &                            (ndiv(m),m=1,3)
	    goto 12
	  endif
	  if(nproc_k>npnt) then
	    if(maswrk) write(iun,*)' NPROC_K > NPNT'
	    call ending
	  endif
c ---------------------------------------------------------------------
	  if(maswrk) then
	    write(iun,'(a7,3i5)')' NDIV: ',(ndiv(m),m=1,3)
	    write(iun,'(a7,3i5)')'NRDIV: ',(nrdiv(m),m=1,3)
        write(iun,'(a9,i5,a9,i5)')'NR_FULL: ',nr_full,' NR_RED: ',nr_red
	  endif
      call timel('*** Divide_BZ   finished ***********')
      end
