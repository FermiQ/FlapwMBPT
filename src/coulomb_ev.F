      subroutine coulomb_ev
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
      integer :: i,j,n_pbt,new,n,iatom,ind,isort
      complex*16 :: c
      real*8, allocatable :: e(:),b(:)
      complex*16, allocatable :: v_tmp(:,:),ovr(:,:),s(:,:),bm(:,:)
      if(allocated(v_opt_e)) deallocate(v_opt_e)
c ------------------------------------------------------------------
      if(me_k==0) then
        n_pbt=n_pbmt+nplwgw(1)
        allocate(v_tmp(n_pbt,n_pbt))
        call v_coul_full(nplwgw(1),v_tmp,1,0)
      endif
      if(me==0) then
        allocate(ovr(n_pbt,n_pbt))
        ovr=(0.d0,0.d0)
        do i=1,n_pbmt
          ovr(i,i)=(1.d0,0.d0)
        enddo
        allocate(s(nplwgw(1),nplwgw(1)))
        call s_overlap(0,s,nplwgw(1),nplwgw(1),indpw_gw(1,1))
        s=s/amega
        do j=1,nplwgw(1)
          do i=1,nplwgw(1)
            ovr(n_pbmt+i,n_pbmt+j)=s(i,j)
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
            ovr(i,j)=ovr(i,j)-conjg(pw_pb(i,1))*pw_pb(j,1)
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
            ind=iopb(iatom)
            c=c+v_tmp(ind,j)*pw_pb(ind,1)
          enddo
          do i=n_pbmt+1,n_pbt
            c=c+v_tmp(i,j)*pw_pb(i,1)
          enddo
          do i=1,n_pbt
            bm(j,i)=conjg(v_tmp(i,j))
          enddo
          do iatom=1,natom
            isort=is(iatom)
            ind=iopb(iatom)
            bm(j,ind)=bm(j,ind)-c*sqrt(4.d0*pi*smt(isort)**3/3.d0/amega)
          enddo
          bm(j,n_pbmt+1)=bm(j,n_pbmt+1)-c
        enddo
        n_opt_pb=new
      endif
#ifdef MPI
      if(nproc/=1) call brdcst(msgint,n_opt_pb,4,master,MPI_COMM_WORLD)
#endif
      allocate(v_opt_e(n_opt_pb))
      allocate(bm_coef(n_opt_pb,n_pbmt+nplwgw(1)))
      if(me==0) then
        v_opt_e=e(1:n_opt_pb)
        call zcopy(n_opt_pb*(n_pbmt+nplwgw(1)),bm,1,bm_coef,1)
        deallocate(e,bm)
      endif
      if(me_k==0) deallocate(v_tmp)
      if(nproc/=1) then
#ifdef MPI
        call brdcst(msgdbl,v_opt_e,8*n_opt_pb,master,MPI_COMM_WORLD)
        call brdcst(msgdbl,bm_coef,16*n_opt_pb*(n_pbmt+nplwgw(1)),
     &              master,MPI_COMM_WORLD)
#endif
      endif
      end
