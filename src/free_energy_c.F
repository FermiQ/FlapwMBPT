      subroutine heg_free_energy_c
c - Evaluates  - Tr ( Ln G  - Ln Gx )   ------------------
      use etot_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: i_omega,i,ispin
	real*8 :: w,et1,d,wg
	complex*16 :: gx,p1,p2,sigc
	d_free_c=0.d0
	do ispin=1,nspin
        do i=0,n_k_heg
          wg=k_m_heg(i)**2*k_w_heg(i)/pi/pi/ro_el_gas/nspin
	    et1=0.d0
	    do i_omega=0,n_omega
	      w=w_omega(i_omega)
	      sigc=sigc_k_omega_heg(i,i_omega,ispin)
	      gx=(1.d0,0.d0)/(dcmplx(0.d0,w)+chem_pot-e_star_x(i,ispin))
	      p1=(1.d0,0.d0)-gx*sigc
	      p2=conjg(p1)
	      d=real(p1*p2)
	      d=log(d)
            et1=et1+et_sum_sc_gx(i_omega)*d
          enddo   !!! over ind_omega
          d_free_c=d_free_c-wg*et1
        enddo
      enddo
      end



      subroutine free_energy_c0(ispin,ind_k,sig,n)
c -   Evaluates  - Tr ( Ln G  - Ln Gx )   ------------------
      use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: ispin,ind_k,n
      real*8, intent(in) :: sig(n,n,2,ndim3_omega)
	  integer :: i_omega,i,ind_omega,k
	  real*8 :: w,et1,sl
	  real*8, allocatable :: d(:)
      complex*16, allocatable :: p1(:,:),p2(:,:),sigc(:,:),g(:,:),gx(:)
	  k=n3_mpi_k(me_k+1)+ind_k
	  allocate(p1(n,n))
	  allocate(p2(n,n))
	  allocate(sigc(n,n))
	  allocate(g(n,n))
	  allocate(d(n))
	  allocate(gx(n))
	  et1=0.d0
	  do ind_omega=1,ndim3_omega
	    i_omega=me_t*ndim3_omega+ind_omega-1
	    w=w_omega(i_omega)
	    call ferm_unpack_omega(g,sig(1,1,1,ind_omega),n,n,n)
	    do i=1,n
	      gx(i)=(1.d0,0.d0)/(dcmplx(0.d0,w)+chem_pot-e_bnd(i,k,ispin))
	    enddo
        call g_c_omega(1,gx,g,n,n)
	    p1=transpose(conjg(g))
        call zgemm('n','n',n,n,n,(1.d0,0.d0),g,n,p1,n,(0.d0,0.d0),p2,n)
	    if(me_b==0) call eig_val_solver(n,n,p2,d)
        if(nproc_b/=1) call brdcst(msgdbl,d,8*n,0,comm_b)
        sl=0.d0
	    do i=1,n
	      sl=sl+log(d(i))-log(gx(i)*conjg(gx(i)))
	    enddo
        et1=et1+et_sum_sc_gx(i_omega)*sl
      enddo   !!! over ind_omega
	  if(nproc_t/=1) call DGOP(et1,1,'  +',comm_t)
      d_free_c=d_free_c+wgt(k)*et1*2.d0/nspin/nrel
	  deallocate(p1,p2,sigc,d,g,gx)
      end
