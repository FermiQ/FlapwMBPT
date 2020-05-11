      subroutine v_coulomb(vmm,vmi,vii,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: khbr
      real*8, intent(out) :: vmm(ncmpl*nd_tb_pbmt(me_t_b+1),
     &                           ndim3_k(me_k+1)),
     &                       vii(nplw_gw,ndim3_k(me_k+1)),
     &                       vmi(ncmpl*ndim_tb_n_pbmt(me_t_b+1),nplw_gw,
     &                           ndim3_k(me_k+1))
      integer :: ind_q,iq,n_pbt,j,i,ind,n0,n1,ij
      real*8 :: pi2a,pi8,q(3),q2
      real*8, allocatable :: v0(:,:)
      complex*16, allocatable :: v_tmp(:,:)
      pi2a=(pi+pi)/par
      pi8=8.d0*pi
      n0=n0_tb_pbmt(me_t_b+1)
      n1=n0+nd_tb_pbmt(me_t_b+1)
      do ind_q=1,ndim3_k(me_k+1)
        iq=n3_mpi_k(me_k+1)+ind_q
        n_pbt=n_pbmt+nplwgw(iq)
c ------------------------------ MM --------------------------------
        allocate(v_tmp(n_pbmt,n_pbmt))
        call v_full_new_mm(v_tmp,iq,khbr)
        if(ncmpl==1) then
          allocate(v0(n_pbmt,n_pbmt))
          call pb_r_from_c_mm(v_tmp,v0,pnt(1,iq))
          ij=0
          do j=1,n_pbmt
            do i=1,j
              ij=ij+1
              if(ij<=n0.or.ij>n1) cycle
              vmm(ij-n0,ind_q)=v0(i,j)
            enddo
          enddo
          deallocate(v0)
        else
          ij=0
          do j=1,n_pbmt
            do i=1,j
              ij=ij+1
              if(ij<=n0.or.ij>n1) cycle
              call zcopy(1,v_tmp(i,j),1,vmm(2*(ij-n0)-1,ind_q),1)
            enddo
          enddo
        endif
        deallocate(v_tmp)
c ------------------------------ MI --------------------------------
        if(ndim_tb_n_pbmt(me_t_b+1)/=0) then
          if(ncmpl==1) then
            call v_coul_mi_r(iq,vmi(1,1,ind_q),nplwgw(iq),khbr)
          else
            call v_coul_mi(iq,vmi(1,1,ind_q),nplwgw(iq),khbr)
          endif
        endif
c ------------------------------ II --------------------------------
        do i=1,nplwgw(iq)
          ind=indpw_gw(i,iq)
          q=(gbs(:,ind)+pnt(:,iq))*pi2a
          q2=dot_product(q,q)
          if(khbr==0.and.q2<1.d-12) then
            vii(i,ind_q)=0.d0
          else
            if(khbr==0) vii(i,ind_q)=pi8/q2
            if(khbr==1) vii(i,ind_q)=pi8/(q2+omega_hse**2)
          endif
        enddo  !! over i
      enddo   !! over ind_q
      call timel('**** V_COULOMB finished ************')
      end
