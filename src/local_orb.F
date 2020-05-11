      subroutine local_orb(work,p0,p1,p2,q0,q1,q2,dr,n,fi,dfi,fidot,
     &	                   dfidot,fidot2,dfidot2,adet,tmp,tmp2,c2,h,
     &                     avv,adv,anorm,fidfid,fid2fid2,dtmp,dtmp2,
     &                     dp0,dp1,dp2,dq0,dq1,dq2)
	  use manager_mod
	  use units_mod
      implicit none
	  integer, intent(in) :: n
	  real*8, intent(in) :: p0(0:n),q0(0:n),p1(0:n),fidfid,fid2fid2,
     &	                    q1(0:n),p2(0:n),q2(0:n),fi,dfi,fidot,
     &	                    dfidot,fidot2,dfidot2,adet,dr(0:n),c2,h,
     &                      dp0(0:n),dq0(0:n),dp1(0:n),
     &                      dq1(0:n),dp2(0:n),dq2(0:n)
      real*8, intent(out) :: tmp(0:n),tmp2(0:n),avv,adv,anorm,dtmp(0:n),
     &                       dtmp2(0:n)
	  integer :: ir
	  real*8 :: fidfid2,dqdall
      real*8, allocatable :: work(:)
      allocate(work(0:n))
c ------- Normalization, zero value and slope on MT ------
      avv=-(fidot2*dfidot-dfidot2*fidot)/adet
      adv=-(dfidot2*fi-fidot2*dfi)/adet
c -------- Normalization --------------------------------
	  do ir=0,n
	    work(ir)=p1(ir)*p2(ir)*dr(ir)
      enddo
	  if(irel.ge.1) then
	    do ir=0,n
	      work(ir)=work(ir)+q1(ir)*q2(ir)*dr(ir)/c2
        enddo
	  endif
      fidfid2=dqdall(h,work,n)
	  anorm=avv*avv+adv*adv*fidfid+fid2fid2+2.d0*adv*fidfid2
	  anorm=1.d0/sqrt(anorm)
      do ir=0,n
        tmp(ir)=anorm*(avv*p0(ir)+adv*p1(ir)+p2(ir))
        dtmp(ir)=anorm*(avv*dp0(ir)+adv*dp1(ir)+dp2(ir))
      enddo
      if(irel.ge.1) then
        do ir=0,n
          tmp2(ir)=anorm*(avv*q0(ir)+adv*q1(ir)+q2(ir))
          dtmp2(ir)=anorm*(avv*dq0(ir)+adv*dq1(ir)+dq2(ir))
        enddo
      endif
      deallocate(work)
      end
