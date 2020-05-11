      subroutine mix_up_down(n,n0,up,upo,dn,dno,t)
	  use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
	  use units_mod
	  use vertex_mod
	  implicit none
	  integer, intent(in) :: n,n0
	  real*8, intent(in) :: upo(n0,n0),dno(n0,n0)
      complex*16, intent(in) :: t(n,n)
	  real*8, intent(inout) :: up(n0,n0),dn(n0,n0)
	  real*8, allocatable :: d(:,:,:),s(:,:,:),v2p(:,:),v2po(:,:)
      complex*16, allocatable :: v2(:,:),v2o(:,:),tmp(:,:)
	  allocate(d(n,n,2))
	  allocate(s(n,n,2))
	  allocate(v2(n,n))
	  allocate(v2o(n,n))
	  allocate(v2p(n,n))
	  allocate(v2po(n,n))
	  allocate(tmp(n,n))
	  call unpack_hermit(v2,dn,n,n0,n,(0.d0,0.d0),(1.d0,0.d0))
	  call zgemm('n','c',n,n,n,(1.d0,0.d0),v2,n,t,n,(0.d0,0.d0),tmp,n)
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2,n)
      call pack_hermit(v2,v2p,n,n,n,0.d0,1.d0)
	  call unpack_hermit(v2o,dno,n,n0,n,(0.d0,0.d0),(1.d0,0.d0))
      call zgemm('n','c',n,n,n,(1.d0,0.d0),v2o,n,t,n,(0.d0,0.d0),tmp,n)
      call zgemm('n','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2o,n)
      call pack_hermit(v2o,v2po,n,n,n,0.d0,1.d0)
	  d(:,:,1)=0.5d0*(up(1:n,1:n)-v2p)
	  d(:,:,2)=0.5d0*(upo(1:n,1:n)-v2po)
	  s(:,:,1)=0.5d0*(up(1:n,1:n)+v2p)
	  s(:,:,2)=0.5d0*(upo(1:n,1:n)+v2po)
      d(:,:,1)=adspin*d(:,:,1)+(1.d0-adspin)*d(:,:,2)
      s(:,:,1)=adm_sig*s(:,:,1)+(1.d0-adm_sig)*s(:,:,2)
	  up(1:n,1:n)=s(:,:,1)+d(:,:,1)
	  v2p=s(:,:,1)-d(:,:,1)
	  call unpack_hermit(v2,v2p,n,n,n,(0.d0,0.d0),(1.d0,0.d0))
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),v2,n,t,n,(0.d0,0.d0),tmp,n)
	  call zgemm('c','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2,n)
      call pack_hermit(v2,dn,n,n,n0,0.d0,1.d0)
      deallocate(d,s,v2,v2o,v2p,v2po,tmp)
      end
