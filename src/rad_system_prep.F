      subroutine rad_system_prep(bt,n,rad,key,mm,em,dalf,alf,l,kap,v)
c	key=1 - Exchange is included
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,key,kap,l
      real*8, intent(in) :: rad(0:n),bt,v(0:n),dalf(0:n),alf(0:n)
      real*8, intent(out) :: mm(n,0:n+2),em(n,0:n+2)
	integer :: i,k
      real*8 :: r2,rll,rb
      if(irel==0) then
        do i=1,n
          rb=rad(i)+bt
          r2=2.d0*bt/rad(i)+1.d0
          rll=v(i)*rb*rb+l*(l+1)*(1.d0+bt/rad(i))**2
          do k=0,n+2
            mm(i,k)=-pp012(i,k,2)-r2*pp012(i,k,1)+rll*pp012(i,k,0)
            if(key==1) mm(i,k)=mm(i,k)+rb*rb*f_xc_spl(k,i)
            em(i,k)=rb*rb*pp012(i,k,0)
          enddo
        enddo
      else if(irel==1) then
        do i=1,n
          rb=rad(i)+bt
          r2=2.d0*bt/rad(i)+1.d0-dalf(i)*rb/alf(i)
          rll=v(i)*rb*rb*alf(i)+l*(l+1)*(1.d0+bt/rad(i))**2
          do k=0,n+2
            mm(i,k)=-pp012(i,k,2)-r2*pp012(i,k,1)+rll*pp012(i,k,0)
            if(key==1) mm(i,k)=mm(i,k)+rb*rb*f_xc_spl(k,i)*alf(i)
            em(i,k)=rb*rb*alf(i)*pp012(i,k,0)
          enddo
        enddo
      else if(irel==2) then
        do i=1,n
          rb=rad(i)+bt
          r2=2.d0*bt/rad(i)+1.d0-dalf(i)*rb/alf(i)
          rll=v(i)*rb*rb*alf(i)+kap*(kap+1)*(1.d0+bt/rad(i))**2
     &       +dalf(i)*(kap+1)*rb*rb/alf(i)/rad(i)
          do k=0,n+2
            mm(i,k)=-pp012(i,k,2)-r2*pp012(i,k,1)+rll*pp012(i,k,0)
            if(key==1) mm(i,k)=mm(i,k)+rb*rb*f_xc_spl(k,i)*alf(i)
            em(i,k)=rb*rb*alf(i)*pp012(i,k,0)
          enddo
        enddo
      endif
      end