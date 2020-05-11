      subroutine integral_pbmt_pw_0_red(mt_pw_0,qq)
	use atom_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: qq
      complex*16, intent(out) :: mt_pw_0(maxlfpb_red,nsort)
	integer :: isort,l,ie,i,n,mtpb,ir
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
        allocate(iik(n,0:3,0:lmpb_red(isort)))
        allocate(f(0:n,2))
        call integral_m_j_prep(n,lmpb_red(isort),gg,r(0,isort),iik)
        do l=0,lmpb_red(isort)
          cc=(0.d0,-1.d0)**l*pi4
	    do ie=1,ntle_pb_red(l,isort)
	      i=ind_prod0_red(ie,l,isort)
	      mtpb=ind_pb_red(i,isort)
            do ir=0,n
              f(ir,1)=fun_pb_red(mtpb+ir)*r(ir,isort)**2
              f(ir,2)=fun_pb_red(mtpb+ir)*2*r(ir,isort)
     &               +dfun_pb_red(mtpb+ir)*r(ir,isort)**2
            enddo
            call integral_m_j_fin(n,f,r(0,isort),iik(1,0,l),res)
            mt_pw_0(i,isort)=cc*res
	    enddo  !! over ie
        enddo !! over l
        deallocate(iik,f)
      enddo   !! over isort
      end