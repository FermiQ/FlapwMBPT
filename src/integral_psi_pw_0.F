      subroutine integral_psi_pw_0(mt_pw_0,qq,ispin)
	  use atom_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      real*8, intent(in) :: qq
      complex*16, intent(out) :: mt_pw_0(maxlfun,nsort)
	  integer :: isort,l,ie,n,ir,it,jn,nres,mtw
      real*8 :: pi4,res,pi2a,gg,sm,b
      complex*16 :: cc
      real*8, allocatable :: iik(:,:,:),f(:,:)
	  pi4=4.d0*pi
      pi2a=(pi+pi)/par
      gg=pi2a*qq
	  do isort=1,nsort
        n=nrad(isort)
        sm=smt(isort)
        b=b_rad(isort)
        allocate(iik(n,0:3,0:lmb(isort)))
        allocate(f(0:n,2))
        call integral_m_j_prep(n,lmb(isort),gg,r(0,isort),iik)
        do l=0,lmb(isort)
          cc=(0.d0,-1.d0)**l*pi4
	      do ie=1,ntle(l,isort)
            it=1
            if(augm(ie,l,isort)=='LAW') it=2
            if(augm(ie,l,isort)=='AWL') it=2
            do jn=1,it
              nres=indfun0(jn,ie,l+1,isort)
              mtw=ind_wf(nres,isort)
              do ir=0,n
                f(ir,1)=gfun(mtw+ir,ispin)*r(ir,isort)**2
                f(ir,2)=gfun(mtw+ir,ispin)*2*r(ir,isort)
     &                 +gfun1(mtw+ir,ispin)*r(ir,isort)**2
              enddo
              call integral_m_j_fin(n,f,r(0,isort),iik(1,0,l),res)
              mt_pw_0(nres,isort)=cc*res
            enddo
	      enddo  !! over ie
        enddo !! over l
        deallocate(iik,f)
      enddo   !! over isort
      end
