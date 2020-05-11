      subroutine eqr_vertex_wan(eqr)
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      complex*16, intent(out) :: eqr(nqdiv_c,ndim_k_red(me_k+1))
	integer :: ir,iq,ind_q
	real*8 :: pi2,fas
	pi2=pi+pi
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
	  do ir=1,nqdiv_c
	    fas=pi2*dot_product(pnt_c(:,iq),rr_coarse(:,ir))
	    eqr(ir,ind_q)=dcmplx(cos(fas),sin(fas))
	  enddo
	enddo
      end
      
      
      
      subroutine eqr_vertex_line(eqr,iq)
	use atom_mod
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: iq
      complex*16, intent(out) :: eqr(nqdiv_c)
	integer :: ir
	real*8 :: pi2,fas
	pi2=pi+pi
	do ir=1,nqdiv_c
	  fas=pi2*dot_product(pnt_c(:,iq),rr_coarse(:,ir))
	  eqr(ir)=dcmplx(cos(fas),sin(fas))
	enddo
      end