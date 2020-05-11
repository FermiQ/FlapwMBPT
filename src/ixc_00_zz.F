      subroutine ixc_00_zz(ind_nu,nn,chi0,ixc0,chi)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,nn
      complex*16, intent(in) :: chi0(nn,nn),ixc0(nn,nn)
      complex*16, intent(out) :: chi(nn,nn)
      integer :: n_pbt,i,nne,n3,info,n6,n36
      real*8 :: pi2a,pi8
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: tmp(:,:),bb(:,:),cc(:,:),rsp(:,:),
     &                           vt(:,:),pt(:,:)
      pi8=8.d0*pi
      pi2a=(pi+pi)/par
      n_pbt=n_pbmt+nplwgw(1)
      nne=n_opt_pb
      n3=nne+3
      n6=n3+n3
      allocate(tmp(n_pbt,n_pbt))
      allocate(bb(n_pbt,n_pbt))
      allocate(cc(n_pbt,n_pbt))
      n_ixc_0=2*nspin-1
      if(nspin==1) then
        allocate(rsp(n3,n3))
c -------- Transform P to the orthonormal basis set ----------------
        bb=chi0
        call zgemm('n','c',n_pbt,nne,n_pbt,(1.d0,0.d0),bb,n_pbt,
     &             bm_coef,nne, (0.d0,0.d0),tmp,n_pbt)
        call zgemm('n','n',nne,nne,n_pbt,(1.d0,0.d0),bm_coef,nne,
     &             tmp,n_pbt,(0.d0,0.d0),rsp(4,4),n3)
        rsp(1:3,1:3)=p_head_dft(:,:,ind_nu)
        rsp(1:3,4:n3)=p_wing_dft(:,:,ind_nu)
        rsp(4:n3,1:3)=conjg(transpose(p_wing_dft(:,:,ind_nu)))
        allocate(vt(n3,n3))
        vt=(0.d0,0.d0)
        do i=1,3
          vt(i,i)=pi8
        enddo
c        call zgemm('n','n',n_pbt,nne,n_pbt,(1.d0,0.d0),ixc0,
c     &             n_pbt,v_opt_pb,n_pbt,(0.d0,0.d0),tmp,n_pbt)
c        call zgemm('c','n',nne,nne,n_pbt,(1.d0,0.d0),v_opt_pb,n_pbt,
c     &             tmp,n_pbt,(0.d0,0.d0),vt(4,4),n3)
        do i=1,nne
          vt(i+3,i+3)=vt(i+3,i+3)+v_opt_e(i)
        enddo
        allocate(pt(n3,n3))
        call zgemm('n','n',n3,n3,n3,(1.d0,0.d0),rsp,n3,vt,n3,
     &             (0.d0,0.d0),pt,n3)
        pt=-pt
        do i=1,n3
          pt(i,i)=(1.d0,0.d0)+pt(i,i)
        enddo
        allocate(ipiv(n3))
        call zgesv(n3,n3,pt,n3,ipiv,rsp,n3,info)
        deallocate(vt,pt,ipiv)
        n36=n3
      else
        allocate(rsp(n6,n6))
c -------- Transform P to the orthonormal basis set ----------------
        bb=chi0(1:n_pbt,1:n_pbt)
        call zgemm('n','c',n_pbt,nne,n_pbt,(1.d0,0.d0),bb,n_pbt,
     &             bm_coef,nne, (0.d0,0.d0),tmp,n_pbt)
        call zgemm('n','n',nne,nne,n_pbt,(1.d0,0.d0),bm_coef,nne,
     &             tmp,n_pbt,(0.d0,0.d0),rsp(4,4),n6)
        rsp(1:3,1:3)=p_head_dft(:,:,ind_nu)
        rsp(1:3,4:n3)=p_wing_dft(:,:,ind_nu)
        rsp(4:n3,1:3)=conjg(transpose(p_wing_dft(:,:,ind_nu)))
        rsp(n3+1:n6,n3+1:n6)=rsp(1:n3,1:n3)
        bb=chi0(1:n_pbt,n_pbt+1:nn)
        call zgemm('n','c',n_pbt,nne,n_pbt,(1.d0,0.d0),bb,n_pbt,
     &             bm_coef,nne, (0.d0,0.d0),tmp,n_pbt)
        call zgemm('n','n',nne,nne,n_pbt,(1.d0,0.d0),bm_coef,nne,
     &             tmp,n_pbt,(0.d0,0.d0),rsp(n3+4,4),n6)
        rsp(1:3,n3+1:n3+3)=p_head_0z(:,:,ind_nu)
        rsp(1:3,n3+4:n6)=p_wing_0z(:,:,ind_nu)
        rsp(4:n3,n3+1:n3+3)=conjg(transpose(p_wing_0z(:,:,ind_nu)))
        rsp(n3+1:n6,1:n3)=rsp(1:n3,n3+1:n6)
        allocate(vt(n6,n6))
        vt=(0.d0,0.d0)
        do i=1,3
          vt(i,i)=pi8
        enddo
        bb=ixc0(1:n_pbt,1:n_pbt)
c        call zgemm('n','n',n_pbt,nne,n_pbt,(1.d0,0.d0),bb,
c     &             n_pbt,v_opt_pb,n_pbt,(0.d0,0.d0),tmp,n_pbt)
c        call zgemm('c','n',nne,nne,n_pbt,(1.d0,0.d0),v_opt_pb,n_pbt,
c     &             tmp,n_pbt,(0.d0,0.d0),vt(4,4),n6)
        do i=1,nne
          vt(i+3,i+3)=vt(i+3,i+3)+v_opt_e(i)
        enddo
        bb=ixc0(1:n_pbt,n_pbt+1:nn)
c        call zgemm('n','n',n_pbt,nne,n_pbt,(1.d0,0.d0),bb,
c     &             n_pbt,v_opt_pb,n_pbt,(0.d0,0.d0),tmp,n_pbt)
c        call zgemm('c','n',nne,nne,n_pbt,(1.d0,0.d0),v_opt_pb,n_pbt,
c     &             tmp,n_pbt,(0.d0,0.d0),vt(4,n3+4),n6)
        vt(n3+1:n6,1:n3)=vt(1:n3,n3+1:n6)
        bb=ixc0(n_pbt+1:nn,n_pbt+1:nn)
c        call zgemm('n','n',n_pbt,nne,n_pbt,(1.d0,0.d0),bb,
c     &             n_pbt,v_opt_pb,n_pbt,(0.d0,0.d0),tmp,n_pbt)
c        call zgemm('c','n',nne,nne,n_pbt,(1.d0,0.d0),v_opt_pb,n_pbt,
c     &             tmp,n_pbt,(0.d0,0.d0),vt(n3+4,n3+4),n6)
        allocate(pt(n6,n6))
        call zgemm('n','n',n6,n6,n6,(1.d0,0.d0),rsp,n6,vt,n6,
     &             (0.d0,0.d0),pt,n6)
        pt=-pt
        do i=1,n6
          pt(i,i)=(1.d0,0.d0)+pt(i,i)
        enddo
        allocate(ipiv(n6))
        call zgesv(n6,n6,pt,n6,ipiv,rsp,n6,info)
        deallocate(vt,pt,ipiv)
        n36=n6
      endif
c      call zgemm('n','c',nne,n_pbt,nne,(1.d0,0.d0),rsp(4,4),n36,
c     &           v_opt_pb,n_pbt,(0.d0,0.d0),bb,n_pbt)
c      call zgemm('n','n',n_pbt,n_pbt,nne,(1.d0,0.d0),v_opt_pb,n_pbt,
c     &           bb,n_pbt,(0.d0,0.d0),tmp,n_pbt)
      chi(1:n_pbt,1:n_pbt)=tmp
      p_head_rsp(:,:,ind_nu)=rsp(1:3,1:3)
      if(nspin==2) then
c        call zgemm('n','c',nne,n_pbt,nne,(1.d0,0.d0),rsp(4,n3+4),n36,
c     &             v_opt_pb,n_pbt,(0.d0,0.d0),bb,n_pbt)
c        call zgemm('n','n',n_pbt,n_pbt,nne,(1.d0,0.d0),v_opt_pb,n_pbt,
c     &             bb,n_pbt,(0.d0,0.d0),tmp,n_pbt)
        chi(1:n_pbt,n_pbt+1:nn)=tmp
        chi(n_pbt+1:nn,1:n_pbt)=tmp
c        call zgemm('n','c',nne,n_pbt,nne,(1.d0,0.d0),rsp(n3+4,n3+4),n36,
c     &             v_opt_pb,n_pbt,(0.d0,0.d0),bb,n_pbt)
c        call zgemm('n','n',n_pbt,n_pbt,nne,(1.d0,0.d0),v_opt_pb,n_pbt,
c     &             bb,n_pbt,(0.d0,0.d0),tmp,n_pbt)
        chi(n_pbt+1:nn,n_pbt+1:nn)=tmp
      endif
      deallocate(rsp,tmp,bb,cc)
      end
