      subroutine dp_v_vrt
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i_nu,j,i,ind_nu,ind_q,iq,nw,np,nne,ij
      real*8 :: de,pi8
      complex*16 :: cc
      complex*16, allocatable :: pc(:,:),p11(:),tmp(:,:)
      pi8=8.d0*pi
      nw=n_pbtot_red
      nne=n_opt_pb_red
      e_vertex_pw_tau0=0.d0
      allocate(pc(nw,nw))
      do ind_q=1,ndim3_k(me_k+1)
        iq=n3_mpi_k(me_k+1)+ind_q
        np=n_pbmt_red+nplwgw_red(iq)
        if(iq==1) allocate(p11(nne))
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          call unpack_hermit(pc,p_vrt_nu_q(1,1,ind_nu,ind_q),np,nw,nw,
     &                       (0.d0,0.d0),(1.d0,0.d0))
          if(iq==1) then
            allocate(tmp(np,nne))
            call zgemm('n','c',np,nne,np,(1.d0,0.d0),pc,nw,bm_coef_red,
     &                 nne,(0.d0,0.d0),tmp,np)
            p11=(0.d0,0.d0)
            do i=1,nne
              do j=1,np
                p11(i)=p11(i)+bm_coef_red(i,j)*tmp(j,i)
              enddo
            enddo
            deallocate(tmp)
            if(metal.and.i_nu==0) then
              call pw_from_pb_red(1,np,nw,pc,cc,3)
              de=pi8*q2aver_c*cc
            else
              cc=(0.d0,0.d0)
              do i=1,3
                ij=ind_tens33(i,i)
                cc=cc+p_head_vrt(ij,ind_nu)
              enddo
              de=pi8*cc/3.d0
            endif
            do i=1,nne
              de=de+p11(i)*v_opt_e_red(i)
            enddo
          else
            de=0.d0
            do j=1,np
              do i=1,np
                de=de+pc(i,j)*v_red_q(j,i,iq)
              enddo
            enddo
          endif
          e_vertex_pw_tau0=e_vertex_pw_tau0
     &                    -et_sum_nu(i_nu)*wgt(iq)*de
        enddo  !! over ind_nu
        if(iq==1) deallocate(p11)
      enddo  !! over ind_q
      deallocate(pc)
      if(nproc_t/=1) call dgop(e_vertex_pw_tau0,1,'  +',comm_t)
      if(nproc_k/=1) call dgop(e_vertex_pw_tau0,1,'  +',comm_k)
      end
