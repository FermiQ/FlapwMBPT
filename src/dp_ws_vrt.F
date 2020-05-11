      subroutine dp_ws_vrt
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i_nu,j,i,ind_nu,ind_q,iq,nw,np,nne,ij
      real*8 :: de,p0
      complex*16 :: cc,p33(3,3),w33(3,3)
      complex*16, allocatable :: pc(:,:),p01(:),ws01(:),p11(:,:),
     &                           tmp(:,:),pwing(:,:),wwing(:,:)
      nw=n_pbtot_red
      nne=n_opt_pb_red
      psi_p_vrt_ws=0.d0
      allocate(pc(nw,nw))
      do ind_q=1,ndim3_k(me_k+1)
        iq=n3_mpi_k(me_k+1)+ind_q
        np=n_pbmt_red+nplwgw_red(iq)
        if(iq==1) allocate(p11(nne,nne))
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          call unpack_hermit(pc,p_vrt_nu_q(1,1,ind_nu,ind_q),np,nw,nw,
     &                       (0.d0,0.d0),(1.d0,0.d0))
          if(iq==1) then
            allocate(tmp(np,nne))
            call zgemm('n','c',np,nne,np,(1.d0,0.d0),pc,nw,bm_coef_red,
     &                 nne,(0.d0,0.d0),tmp,np)
            call zgemm('n','n',nne,nne,np,(1.d0,0.d0),bm_coef_red,nne,
     &                 tmp,np,(0.d0,0.d0),p11,nne)
            deallocate(tmp)
            if(metal.and.i_nu==0) then
              allocate(p01(nne))
              p01=(0.d0,0.d0)
              do j=1,np
                cc=(0.d0,0.d0)
                do i=1,np
                  cc=cc+conjg(pw_pb_red(i,1))*pc(i,j)
                enddo
                do i=1,nne
                  p01(i)=p01(i)+cc*conjg(bm_coef_red(i,j))
                enddo
              enddo
              allocate(ws01(nne))
              ws01=(0.d0,0.d0)
              if(ncmpl==1) then
                do j=1,nne
                  do i=1,n_opt_pb
                    ws01(j)=ws01(j)+ws01_met_r(i)*bb_red(i,j)
                  enddo
                enddo
              else
                do j=1,nne
                  do i=1,n_opt_pb
                    ws01(j)=ws01(j)+ws01_met(i)*bb_red(i,j)
                  enddo
                enddo
              endif
              call pw_from_pb_red(1,np,nw,pc,cc,3)
              p0=cc
              call p_ws_conv_q0_met(nne,nne,nne,p01,ws01,p0,
     &                              ws_head_nu(0),p11,
     &                              w11_red(1,1,ind_nu),de,q2aver)
              deallocate(p01,ws01)
            else
              allocate(pwing(3,nne))
              allocate(wwing(3,nne))
              pwing=(0.d0,0.d0)
              wwing=(0.d0,0.d0)
              do j=1,3
                do i=1,3
                  ij=ind_tens33(i,j)
                  p33(i,j)=p_head_vrt(ij,ind_nu)
                enddo
              enddo
              if(ncmpl==1) then
                w33=ws_head_tens_r(:,:,ind_nu)
                call p_ws_conv_q0(nne,nne,nne,pwing,wwing,p33,w33,p11,
     &                            w11_red(1,1,ind_nu),de)
              else
                call p_ws_conv_q0(nne,nne,nne,pwing,wwing,p33,
     &                            ws_head_tens(1,1,ind_nu),p11,
     &                            w11_red(1,1,ind_nu),de)
              endif
              deallocate(pwing,wwing)
            endif
          else
            de=0.d0
            do j=1,np
              do i=1,np
                de=de+pc(i,j)*w_red_q(j,i,ind_nu,iq)
              enddo
            enddo
          endif
          psi_p_vrt_ws=psi_p_vrt_ws+et_sum_nu(i_nu)*wgt(iq)*de
        enddo  !! over ind_nu
        if(iq==1) deallocate(p11)
      enddo  !! over ind_q
      deallocate(pc)
      if(nproc_t/=1) call dgop(psi_p_vrt_ws,1,'  +',comm_t)
      if(nproc_k/=1) call dgop(psi_p_vrt_ws,1,'  +',comm_k)
      end
