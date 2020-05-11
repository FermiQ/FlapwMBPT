      subroutine vertex_p_psi(key)
c     key = 0 - Based on DFT or HF Green's function
c     key = 1 - Based on QP Green's function
c     key = 2 - Based on general Green's function
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      character*3 :: appr
      integer :: iter_w,iter_t,ispin,mx2,nt
      real*8, allocatable :: g_rr_tau(:,:,:,:,:,:),wst4loc(:,:,:,:,:)
      complex*16, allocatable :: tmatr(:,:,:,:,:,:)
      if(psi_p=='00000') return
      mx2=maxel_red**2
      iter_w=0
      iter_t=0
      if(psi_p(1:1)/='0') iter_w=1
      if(psi_p(3:3)/='0') iter_w=2
      if(psi_p(2:2)/='0') iter_t=1
      if(psi_p(1:1)=='6') appr='bnd'
      if(psi_p(1:1)=='4') appr='wmt'
      if(psi_p(1:1)=='3') appr='vmt'
      if(psi_p(1:1)=='2') appr='loc'
      if(appr=='loc') then
        call vertex_p_loc(iter_w,iter_t)
      else if(psi_p(1:3)/='000') then
        call vertex_p_uni_test(iter_w,iter_t)
        if(appr=='vmt') then
          call vertex_p_vmt(iter_w,iter_t,key)
        else
          call vertex_bnd_gw(iter_w,iter_t,appr,key)
        endif
      endif
      if(psi_p(5:5)/='0') then
        if(psi_p(5:5)=='3') then
          allocate(g_rr_tau(nfun_red,nfun_red,nqdiv_c,2,ndim3_tau,
     &                      nspin))
          do ispin=1,nspin
            call g_real_space_mm(g_rr_tau(1,1,1,1,1,ispin),ispin)
          enddo
          allocate(wst4loc(mx2,mx2,natom,nspin,nspin))
          call w4_loc_stat(wst4loc)
          nt=2
          if(psi_p(2:2)/='0') nt=3
          allocate(tmatr(n_lem2_red,n_lem2_red,ndim3_nu,
     &                   ndimc_kk(me_k+1),nspin,nspin))
          call t_matr_sf_gen(tmatr,g_rr_tau,wst4loc,'ph',nt)
          call p_psi5_g2_w1_stat(tmatr,'ph')
          call t_matr_sf_gen(tmatr,g_rr_tau,wst4loc,'pp',nt)
          call p_psi5_g2_w1_stat(tmatr,'pp')
          deallocate(g_rr_tau,wst4loc,tmatr)
        endif
      endif
      end
