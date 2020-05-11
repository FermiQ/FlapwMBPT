      subroutine heg_p_uniform(spin_flips,dp)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      logical, intent(in) :: spin_flips
      real*8, intent(out) :: dp
      integer :: it,ispin
      real*8 :: conv_1
      real*8, allocatable :: conv_table(:),lambda_stat(:,:),
     &                       k_pw(:,:,:,:)
      complex*16, allocatable :: lambda_dyn(:,:,:)
c ----------------------------------------------------------------
      allocate(conv_table(2:iter_ladder_xi))
	allocate(lambda_stat(nrr_div_c,nspin))
	allocate(k_pw(nrr_div_c,2,ndim3_tau,nspin))
      dp=0.d0
      if(iter_ladder_xi>1) conv_table=0.d0
	allocate(lambda_dyn(nrr_div_c,0:n_omega,nspin))
c ----------------------------------------------------------------
	call vertex_heg_uniform(iter_ladder_xi,lambda_dyn,lambda_stat,
     &                        conv_table,conv_1,spin_flips,
     &                        iter_ladder_xi)
c -------------- Correction to the Polarizability ----------------------
      do ispin=1,nspin
        call k_heg_uniform(ispin,lambda_dyn(1,0,ispin),
     &                     k_pw(1,1,1,ispin),lambda_stat(1,ispin))
	  if(me_t==0) dp=dp-real(k_pw(1,1,1,ispin))
	enddo
	if(nspin==1) dp=2.d0*dp
	deallocate(lambda_dyn)
	if(nproc_t/=1) call dgop(dp,1,'  +',comm_t)
	if(maswrk.and.iter_ladder_xi>1) then
	  write(iun,*)' Convergency of the ladder sequence, Q =',1
	  write(iun,'(a5,8(i3,6x))')'It    Convergency'
	  do it=2,iter_ladder_xi
	    write(iun,'(i3,8e9.2)')it,conv_table(it)
        enddo
      endif
c ------------------------------------------------------------------
      deallocate(k_pw,lambda_stat,conv_table)
      call timel('**** HEG_P_UNIFORM finished ********')
      end