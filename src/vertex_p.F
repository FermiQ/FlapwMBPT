      subroutine vertex_p(w_old,key)
c     key = 0 - Based on DFT or HF Green's function
c     key = 1 - Based on QP Green's function
c     key = 2 - Based on general Green's function
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical, intent(in) :: w_old
      integer, intent(in) :: key
      integer :: ind_q,iq,k
      real*8, allocatable :: p_line(:,:,:,:),ph(:,:),ph0(:,:)
      if(w_old) return
      if(ubi=='psi'.and.psi_p=='00000') return
      if(ubi=='bsp'.and.bse_kernel_p=='0000') return
      if(ubi==' qp'.and.ubi_0==' gw'.and.bse_kernel_p=='0000') return
      if(allocated(p_vertex_nu_q)) deallocate(p_vertex_nu_q)
      allocate(p_vertex_nu_q(n_pbtot_red,n_pbtot_red,ndim3_nu,
     &                       ndim_k_red(me_k+1)))
      p_vertex_nu_q=0.d0
      allocate(w_red_q(n_pbtot_red,n_pbtot_red,ndim3_nu,npnt))
      call w_c_from_b
      allocate(vrt_wi(nbndf_bnd,nbndf_bnd,npnt_c,0:1,nspin))
      vrt_wi=(0.d0,0.d0)
      if(ubi==' qp') call vertex_p_bse(key)
      if(ubi=='psi') call vertex_p_psi(key)
      if(ubi=='bsp') call vertex_p_bse(key)
c ------------------------------------------------------------------
      call p_red_q_interp
      deallocate(p_vertex_nu_q)
c -----------------------------------------------------------------
      allocate(p_line(n_pbtot_red,n_pbtot_red,ndim3_nu,n_line))
      p_line=0.d0
      do ind_q=1,ndim3_k(me_k+1)
        iq=n3_mpi_k(me_k+1)+ind_q
        do k=1,n_line
          if(kline_in_npnt(k)==iq) p_line(:,:,:,k)=
     &                               p_vrt_nu_q(:,:,:,ind_q)
        enddo
      enddo
      if(nproc_k/=1) call dgop(p_line,n_pbtot_red**2*ndim3_nu*n_line,
     &                         '  +',comm_k)
      call output_p_vertex_red('P_VERTEX',p_line)
      deallocate(p_line)
c ------------------------------------------------------------------
      allocate(ph(n_k_head_0,ndim3_nu))
      call p_head_vertex_prep(ph)
      if(.not.allocated(p_head_vrt)) allocate(p_head_vrt(6,ndim3_nu))
      allocate(ph0(n_phead_dir,ndim3_nu))
      call p_head_num(ph,p_head_vrt,ph0)
      deallocate(ph,ph0)
      call dp_v_vrt   ! calculates e_vertex_pw_tau0
      call dp_ws_vrt  ! evaluates psi_p_vrt_ws
      deallocate(w_red_q)
      k=4*nbndf_bnd**2*npnt_c*nspin
      if(nproc_t/=1) call dgop(vrt_wi,k,'  +',comm_t)
      if(nproc_k/=1) call dgop(vrt_wi,k,'  +',comm_k)
      if(nproc_b/=1) call dgop(vrt_wi,k,'  +',comm_b)
      call output_vrt_wi
      deallocate(vrt_wi)
      call timel('**** VERTEX_P finished *************')
      end
