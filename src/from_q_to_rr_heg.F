	subroutine from_q_to_rr_heg_red(n,q,fr)
c     Usage:
c     n = nstar  : for Q_irr ----- >  RR
c     n = nq_sph : for Q_sph ----- >  RR
	use heg_mod
	use manager_mod
	use solid_mod
      use vertex_mod
	implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: q(n)
	real*8, intent(out) :: fr(nrr_div_c)
	integer :: i,iq
	complex*16, allocatable :: cc(:)
	allocate(cc(nrr_div_c))
	cc=(0.d0,0.d0)
	do i=1,nq_sph_c
	  iq=q_sph_in_cube_c(i)
	  if(n==nq_sph_c) cc(iq)=q(i)
	  if(n==nstar_c) cc(iq)=q(i_kref_c(i))
	enddo
	call fft3(nr_div_c(1),nr_div_c(2),nr_div_c(3),1,cc,1)
	fr=cc/nqdiv_c/amega
	deallocate(cc)
	end
	
	
	
	subroutine from_q_to_rr_heg_cube_red(fc)
	use heg_mod
	use manager_mod
	use solid_mod
      use vertex_mod
	implicit none
	complex*16, intent(inout) :: fc(nrr_div_c)
	call fft3(nr_div_c(1),nr_div_c(2),nr_div_c(3),1,fc,1)
	fc=fc/nqdiv_c/amega
	end
	
	
	subroutine from_k_to_r_heg1(fk,fr)
	use heg_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	real*8, intent(in) :: fk(0:n_k_heg)
	real*8, intent(out) :: fr(0:n_r_heg)
	call dgemv('n',n_r_heg+1,n_k_heg+1,1.d0,k_to_r,n_r_heg+1,fk,1,
     &	         0.d0,fr,1)
	end	
	
	subroutine from_q_to_r_heg1(fq,fr)
	use heg_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	real*8, intent(in) :: fq(0:n_q_heg)
	real*8, intent(out) :: fr(0:n_r_heg)
	call dgemv('n',n_r_heg+1,n_q_heg+1,1.d0,q_to_r,n_r_heg+1,fq,1,
     &	         0.d0,fr,1)
	end
