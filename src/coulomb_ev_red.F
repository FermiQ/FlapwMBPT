      subroutine coulomb_ev_red
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer :: i,j,n_pbt,new,n,np,iatom,ind,isort
      complex*16 :: c
      real*8, allocatable :: e(:),b(:)
      complex*16, allocatable :: v_tmp(:,:),ovr(:,:),s(:,:),bm(:,:),
     &                           v_opt_pb_red(:,:)
      if(allocated(v_opt_e_red)) deallocate(v_opt_e_red)
c ------------------------------------------------------------------
      if(me_k==0) then
        n_pbt=n_pbmt_red+nplwgw_red(1)
        allocate(v_tmp(n_pbt,n_pbt))
        v_tmp=v_red_q(1:n_pbt,1:n_pbt,1)
      endif
      if(me==0) then
        allocate(ovr(n_pbt,n_pbt))
        ovr=(0.d0,0.d0)
        do i=1,n_pbmt_red
          ovr(i,i)=(1.d0,0.d0)
        enddo
        allocate(s(nplwgw_red(1),nplwgw_red(1)))
        call s_overlap(0,s,nplwgw_red(1),nplwgw_red(1),
     &                 indpw_gw_red(1,1))
        s=s/amega
        do j=1,nplwgw_red(1)
          do i=1,nplwgw_red(1)
            ovr(n_pbmt_red+i,n_pbmt_red+j)=s(i,j)
          enddo
        enddo
        deallocate(s)
        allocate(e(n_pbt))
        allocate(b(n_pbt))
        allocate(s(n_pbt,n_pbt))
        call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),ovr,n_pbt,
     &             v_tmp,n_pbt,(0.d0,0.d0),s,n_pbt)
        call zgemm('n','c',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),s,n_pbt,ovr,
     &             n_pbt,(0.d0,0.d0),v_tmp,n_pbt)
        deallocate(s)
        do j=1,n_pbt
          do i=1,n_pbt
            ovr(i,j)=ovr(i,j)-conjg(pw_pb_red(i,1))*pw_pb_red(j,1)
          enddo
        enddo
        call eig_val_gener_solver(1.d-5,n_pbt,n_pbt,v_tmp,ovr,new,e,b)
        n=0
        do i=1,new
          if(e(i)/e(new)>eps_coul) n=n+1
        enddo
        n=new-n
        new=new-n
        do i=1,new
          e(i)=e(i+n)
          v_tmp(:,i)=v_tmp(:,i+n)
        enddo
        deallocate(ovr,b)
        allocate(bm(new,n_pbt))
        do j=1,new
          c=(0.d0,0.d0)
          do iatom=1,natom
            isort=is(iatom)
            ind=iopb_red(iatom)
            c=c+v_tmp(ind,j)*pw_pb_red(ind,1)
          enddo
          do i=n_pbmt_red+1,n_pbt
            c=c+v_tmp(i,j)*pw_pb_red(i,1)
          enddo
          do i=1,n_pbt
            bm(j,i)=conjg(v_tmp(i,j))
          enddo
          do iatom=1,natom
            isort=is(iatom)
            ind=iopb_red(iatom)
            bm(j,ind)=bm(j,ind)-c*sqrt(4.d0*pi*smt(isort)**3/3.d0/amega)
          enddo
          bm(j,n_pbmt+1)=bm(j,n_pbmt+1)-c
        enddo
        n_opt_pb_red=new
      endif
#ifdef MPI
      if(nproc/=1) call brdcst(msgint,n_opt_pb_red,4,master,
     &                         MPI_COMM_WORLD)
#endif
      allocate(v_opt_e_red(n_opt_pb_red))
      allocate(v_opt_pb_red(n_pbmt_red+nplwgw_red(1),n_opt_pb_red))
      allocate(bm_coef_red(n_opt_pb_red,n_pbmt_red+nplwgw_red(1)))
      if(me==0) then
        v_opt_e_red=e(1:n_opt_pb_red)
        v_opt_pb_red=v_tmp(:,1:n_opt_pb_red)
        bm_coef_red=bm
        deallocate(e,bm)
      endif
      if(me_k==0) deallocate(v_tmp)
      if(nproc/=1) then
#ifdef MPI
        call brdcst(msgdbl,v_opt_e_red,8*n_opt_pb_red,master,
     &              MPI_COMM_WORLD)
        call brdcst(msgdbl,v_opt_pb_red,
     &              16*n_opt_pb_red*(n_pbmt_red+nplwgw_red(1)),
     &              master,MPI_COMM_WORLD)
        call brdcst(msgdbl,bm_coef_red,
     &              16*n_opt_pb_red*(n_pbmt_red+nplwgw_red(1)),
     &              master,MPI_COMM_WORLD)
#endif
      endif
      allocate(bb_red(n_opt_pb,n_opt_pb_red))
      n_pbt=n_pbmt+nplwgw(1)
      np=n_pbmt_red+nplwgw_red(1)
      allocate(bm(n_opt_pb,np))
      if(ncmpl==1) then
        allocate(bm_coef(n_opt_pb,n_pbmt+nplwgw(1)))
        call bm_c_from_r(nplwgw(1),bm_coef,bm_coef_r)
      endif
      call rpb_from_fpb_right(bm,bm_coef,n_opt_pb,n_pbt,np,1)
      if(ncmpl==1) deallocate(bm_coef)
      call zgemm('n','n',n_opt_pb,n_opt_pb_red,np,(1.d0,0.d0),bm,
     &           n_opt_pb,v_opt_pb_red,np,(0.d0,0.d0),bb_red,
     &           n_opt_pb)
      deallocate(bm,v_opt_pb_red)
      end
