      subroutine p_head_vertex_prep(ph)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(out) :: ph(n_k_head_0,ndim3_nu)
      integer :: i,ind_nu,i_nu,ii,ind_q,iq,np
      complex*16 :: cc
      complex*16, allocatable :: work(:,:)
      ph=(0.d0,0.d0)
      do ind_q=1,ndim3_k(me_k+1)
        iq=n3_mpi_k(me_k+1)+ind_q
        ii=0
        do i=1,n_k_head_0
          if(iq==khead0(i)) ii=i
        enddo
        if(ii==0) cycle
        np=n_pbmt_red+nplwgw_red(iq)
        allocate(work(np,np))
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          call unpack_hermit(work,p_vrt_nu_q(1,1,ind_nu,ind_q),np,
     &                       n_pbtot_red,np,(0.d0,0.d0),(1.d0,0.d0))
          call pw_from_pb_red(iq,np,np,work,cc,3)
          ph(ii,ind_nu)=cc
        enddo
        deallocate(work)
      enddo
      if(nproc_k/=1) call dgop(ph,n_k_head_0*ndim3_nu,'  +',comm_k)
      end
