      subroutine heg_p_line_stat(spin_flips,dp)
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
      real*8, intent(out) :: dp(n_line)
      integer :: it,iq,k,ispin
      real*8 :: conv_1
      real*8, allocatable :: conv_table(:,:)
      complex*16, allocatable :: lambda_stat(:,:),k_pw(:,:,:,:),
     &                           lambda_dyn(:,:,:)
c ----------------------------------------------------------------
      allocate(conv_table(2:iter_ladder_xi,3))
	allocate(lambda_stat(nrr_div,nspin))
	allocate(k_pw(nrr_div,2,ndim3_tau,nspin))
      dp=0.d0
      if(iter_ladder_xi>1) conv_table=0.d0
	allocate(lambda_dyn(nrr_div,0:n_omega,nspin))
c ----------------------------------------------------------------
	do iq=1,n_line
	  call vertex_heg_nu0(iter_ladder_xi,iq,lambda_dyn,lambda_stat,
     &                      conv_table,conv_1,spin_flips,iter_ladder_xi,
     &                      q_line(1,iq))
c -------------- Correction to the Polarizability ----------------------
        do ispin=1,nspin
          call k_heg_stat(ispin,q_line(1,iq),lambda_dyn(1,0,ispin),
     &                    k_pw(1,1,1,ispin),lambda_stat(1,ispin))
	    if(me_t==0) dp(iq)=dp(iq)-real(k_pw(1,1,1,ispin))
	  enddo
      enddo   !! over iq
	deallocate(lambda_dyn)
	if(nproc_t/=1) call dgop(dp,n_line,'  +',comm_t)
	if(nspin==1) dp=2.d0*dp
	if(maswrk.and.iter_ladder_xi>1) then
	  do k=1,3
	    write(iun,*)' Convergency of the ladder sequence, Q =',k
	    write(iun,'(a5,8(i3,6x))')'It\Nu',0,1,2,3,4,5,6,7
	    do it=2,iter_ladder_xi
	      write(iun,'(i3,8e9.2)')it,conv_table(it,k)
          enddo
        enddo
      endif
c ------------------------------------------------------------------
      deallocate(k_pw,lambda_stat,conv_table)
      call timel('**** HEG_P_LINE_STAT finished ******')
      end