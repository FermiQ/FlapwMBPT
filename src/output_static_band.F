      subroutine output_static_band(xx,key)
c     key = 0 -   V_H
c     key = 1 -   V_XC
c     key = 2 -   SIGMA_X
c     key = 3 -   SIGMA_QP_EFF
c     key = 4 -   delta V_H_SC
c     key = 5 -   delta SIGMA_X_SC
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: key
      complex*16, intent(in) :: xx(nd_tb_ij(me_t_b+1),ndim3_k(me_k+1),
     &                             nspin)
      character*4 ch2
      character*17 ch1
      integer :: k,i_len,ispin,ind_k,i,ib,n0,nmin,nnn,ks,i0,i1,ii
      real*8 :: ff
      real*8, allocatable :: z_out(:,:)
      nnn=8
      ff=evolt/2.d0
      if(ubi=='dft') ch2='.dft'
      if(ubi==' hf') ch2='_.hf'
      if(ubi==' qp') ch2='_.qp'
      if(ubi==' gw') ch2='_.gw'
      if(ubi=='psi') ch2='.psi'
      if(ubi=='bsp') ch2='.bsp'
      if(key==0) ch1='_V_Hartree_band_k'
      if(key==1) ch1='______V_XC_band_k'
      if(key==2) ch1='___Sigma_X_band_k'
      if(key==3) ch1='_Sg_QP_Eff_band_k'
      if(key==4) ch1='_d_V_Hartr_band_k'
      if(key==5) ch1='_d_Sigma_X_band_k'
      n0=nbndf
      nmin=minval(n_bnd)
      allocate(z_out(nmin,n_sym))
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//ch1//ch2)
      endif
      i0=n0_tb_ij(me_t_b+1)
      i1=i0+nd_tb_ij(me_t_b+1)
      do ispin=1,nspin
        if(maswrk) write(3,*)' Ispin = ',ispin
        z_out=0.d0
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          ks=k_npnt_in_sym(k)
          if(ks==0) cycle
          do i=1,nmin
            ii=i*(i+1)/2
            if(ii<=i0) cycle
            if(ii>i1) cycle
            z_out(i,ks)=xx(ii-i0,ind_k,ispin)
          enddo
        enddo
c --------- Collect from all processes  ----------------------------
        if(nproc_k/=1) call dgop(z_out,nmin*n_sym,'  +',comm_k)
        if(nproc_t*nproc_b/=1) call dgop(z_out,nmin*n_sym,'  +',
     &                                   comm_t_b)
        if(maswrk) then
          write(3,*)' N_SYM K-points from ',1,' to ',min(nnn,n_sym)
          write(3,'(4x,8(4x,a1,5x))')(sym0_label(k),k=1,min(n_sym,nnn))
          do ib=1,nmin
            write(3,'(i4,8(1x,f9.4))')ib,(z_out(ib,k)*ff,k=1,
     &	                                  min(nnn,n_sym))
          enddo
        endif
      enddo
      if(maswrk) close(3) 
      deallocate(z_out)
      end
