      subroutine mix_up_down_st(n,n0,up,upo,dn,dno,t)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n0
      complex*16, intent(in) :: upo(nd_tb_ij(me_t_b+1)),
     &                          dno(nd_tb_ij(me_t_b+1)),t(n,n)
      complex*16, intent(inout) :: up(nd_tb_ij(me_t_b+1)),
     &                             dn(nd_tb_ij(me_t_b+1))
      complex*16, allocatable :: v2(:,:),v2o(:,:),tmp(:,:),v2p(:),
     &                           v2po(:),d(:,:),s(:,:)
      allocate(d(nd_tb_ij(me_t_b+1),2))
      allocate(s(nd_tb_ij(me_t_b+1),2))
      allocate(v2(n,n))
      allocate(v2o(n,n))
      allocate(v2p(nd_tb_ij(me_t_b+1)))
      allocate(v2po(nd_tb_ij(me_t_b+1)))
      allocate(tmp(n,n))
      call unpack_hrm_tb(v2,dn,n,n)
      call zgemm('n','c',n,n,n,(1.d0,0.d0),v2,n,t,n,(0.d0,0.d0),tmp,n)
      call zgemm('n','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2,n)
      call pack_hrm_tb(v2,v2p,n,n,0)
      call unpack_hrm_tb(v2o,dno,n,n)
      call zgemm('n','c',n,n,n,(1.d0,0.d0),v2o,n,t,n,(0.d0,0.d0),tmp,n)
      call zgemm('n','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2o,n)
      call pack_hrm_tb(v2o,v2po,n,n,0)
      d(:,1)=0.5d0*(up-v2p)
      d(:,2)=0.5d0*(upo-v2po)
      s(:,1)=0.5d0*(up+v2p)
      s(:,2)=0.5d0*(upo+v2po)
      d(:,1)=adspin*d(:,1)+(1.d0-adspin)*d(:,2)
      s(:,1)=adm_sig*s(:,1)+(1.d0-adm_sig)*s(:,2)
      up=s(:,1)+d(:,1)
      v2p=s(:,1)-d(:,1)
      call unpack_hrm_tb(v2,v2p,n,n)
      call zgemm('n','n',n,n,n,(1.d0,0.d0),v2,n,t,n,(0.d0,0.d0),tmp,n)
      call zgemm('c','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2,n)
      call pack_hrm_tb(v2,dn,n,n,0)
      deallocate(d,s,v2,v2o,v2p,v2po,tmp)
      end
