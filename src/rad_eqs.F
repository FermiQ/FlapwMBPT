      subroutine rad_eqs
	  use atom_mod
	  use etot_mod
	  use manager_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer :: isort,ispin,mt,mt0,nlmb,it0
	  real*8 :: sqpi4,c2
	  integer, allocatable :: nds(:,:,:)
	  real*8, allocatable :: acc(:,:,:),v0(:),work(:),dff(:,:,:)
      c2=clight**2
      ffsmt=0.d0
	  ffhmt=0.d0
	  sqpi4=sqrt(4.d0*pi)
	  if(maxl_cor>=0) allocate(t_x(nrel*maxb+1,nrel*maxl_cor+1,
     & 	                           0:maxb+maxl_cor))
      allocate(t1_x(0:maxb+maxb,maxlfun,maxlfun,nrel*maxb+1))
	  gfun_old=gfun
	  if(irel>=1) gfund_old=gfund
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	  nlmb=nrel*maxb+1
 	  allocate(acc(maxntle,nlmb,nspin))
 	  allocate(dff(maxntle,nlmb,nspin))
 	  allocate(nds(maxntle,nlmb,nspin))
	  allocate(v0(0:maxnrad))
	  allocate(work(0:maxnrad))
	  do isort=1,nsort
        mt0=indmt(1,isort,1)
        do ispin=1,nspin
          mt=indmt(1,isort,ispin)
          v0(0:nrad(isort))=v_mt_h(mt0:mt0+nrad(isort))
          v0(0:nrad(isort))=v0(0:nrad(isort))
     &                     +v_mt_xc(mt:mt+nrad(isort))
	      v0(0:nrad(isort))=v0(0:nrad(isort))/sqpi4
	      call radli(v0(0),isort,ispin,acc,dff,nds)
	    enddo
	    call prnatom(1,isort,nspin,acc,dff,nds)
	  enddo
	  deallocate(acc,v0,work,dff,nds)
	  if(maxl_cor>=0) deallocate(t_x)
	  deallocate(t1_x)
	  it0=iter_dft
	  if(iter==1.or.iter==it0) then
	    gfun_old=gfun
	    if(irel>=1) gfund_old=gfund
	  endif
      call timel('******* RAD_EQS finished ***********')
      end
