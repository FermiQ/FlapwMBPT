      subroutine integral_psi_pw(mt_pw,mt_pw_0,v)
	use atom_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: v(3)
      complex*16, intent(in) :: mt_pw_0(maxlfun,nsort)
      complex*16, intent(out) :: mt_pw(nfun)
	integer :: iatom,isort,lm,km,ind0,i
      real*8 :: pi4,sqpi4,qq,pi2,phase
      complex*16 :: cc
	real*8, allocatable :: ylm(:)
	allocate (ylm(limlb))
	pi2=pi+pi
	pi4=4.d0*pi
	sqpi4=sqrt(pi4)
	qq=sqrt(dot_product(v,v))
	if(qq.gt.1.d-9) then
	  call sphharm(v(1)/qq,v(2)/qq,v(3)/qq,maxb,ylm)
	else
	  ylm=0.d0
	  ylm(1)=1.d0/sqpi4
	endif
	do iatom=1,natom
	  isort=is(iatom)
	  ind0=io_lem(iatom)-1
	  phase=-pi2*dot_product(v,tau(:,iatom))
	  cc=dcmplx(cos(phase),sin(phase))
        do km=1,lfunm(isort)
	    lm=lm_isz(km,isort)
	    i=lf_isz(km,isort)
          mt_pw(ind0+km)=cc*ylm(lm)*mt_pw_0(i,isort)
	  enddo !! over km
	enddo   !! over iatom
	deallocate (ylm)
      end
      
      
      
      
      subroutine integral_phi_pw(mt_red,mt_pw_0,v)
	use atom_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: v(3)
      complex*16, intent(in) :: mt_pw_0(maxlfun,nsort)
      complex*16, intent(out) :: mt_red(nfun_red)
	integer :: iatom,isort,lm,km,ind0,i,km0
      real*8 :: pi4,sqpi4,qq,pi2,phase
      complex*16 :: cc
	real*8, allocatable :: ylm(:)
	allocate (ylm(limlb))
	pi2=pi+pi
	pi4=4.d0*pi
	sqpi4=sqrt(pi4)
	qq=sqrt(dot_product(v,v))
	if(qq.gt.1.d-9) then
	  call sphharm(v(1)/qq,v(2)/qq,v(3)/qq,maxb,ylm)
	else
	  ylm=0.d0
	  ylm(1)=1.d0/sqpi4
	endif
	do iatom=1,natom
	  isort=is(iatom)
	  ind0=io_lem_red(iatom)-1
	  phase=-pi2*dot_product(v,tau(:,iatom))
	  cc=dcmplx(cos(phase),sin(phase))
        do km0=1,lfunm_red(isort)
          km=lme_red(km0,isort)
	    lm=lm_isz(km,isort)
	    i=lf_isz(km,isort)
          mt_red(ind0+km0)=cc*ylm(lm)*mt_pw_0(i,isort)
	  enddo !! over km
	enddo   !! over iatom
	deallocate (ylm)
      end