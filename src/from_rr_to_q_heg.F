	subroutine from_rr_to_q_heg(q,fr)
	use heg_mod
	use manager_mod
	use solid_mod
	implicit none
	real*8, intent(in) :: fr(nrr_div_c)
	real*8, intent(out) :: q(nstar_c)
	integer :: i,iq,istar
	real*8 :: a
	complex*16, allocatable :: cc(:)
	allocate(cc(nrr_div_c))
	cc=fr
	call fft3(nr_div_c(1),nr_div_c(2),nr_div_c(3),1,cc,-1)
	do istar=1,nstar_c
	  a=0.d0
	  do i=indstar_c(istar-1)+1,indstar_c(istar)
	    iq=q_sph_in_cube_c(i)
	    a=a+cc(iq)
	  enddo
	  q(istar)=a/(indstar_c(istar)-indstar_c(istar-1))
	enddo
	q=q*amega/nr_full_c
	deallocate(cc)
	end
	
	
	
	subroutine from_rr_to_q_heg_cube(fc)
	use heg_mod
	use manager_mod
	use solid_mod
	implicit none
	complex*16, intent(inout) :: fc(nrr_div)
	call fft3(nr_div(1),nr_div(2),nr_div(3),1,fc,-1)
	fc=fc*amega/nr_full
	end
	
	subroutine from_rr_to_q_heg_cube_red(fc)
	use heg_mod
	use manager_mod
	use solid_mod
	implicit none
	complex*16, intent(inout) :: fc(nrr_div_c)
	call fft3(nr_div_c(1),nr_div_c(2),nr_div_c(3),1,fc,-1)
	fc=fc*amega/nr_full_c
	end