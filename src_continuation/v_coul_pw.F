      subroutine v_coul_pw(vii,n1,n2,iq)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n1,n2,iq
      complex*16, intent(out) :: vii
      integer :: j,i,iatom,isort,n
      real*8 :: a,v(3),pi2a,v1,pi2,phase,v0(3)
      complex*16 :: cc
      complex*16, allocatable :: vc1(:),vc2(:)
      pi2=pi+pi
      pi2a=pi2/par
      n=nplwro
      allocate(vc1(n),vc2(n))
      vc1=(0.d0,0.d0)
      i=indpw_gw(n1,iq)
      v0=gbs(:,i)*pi2a
      do j=1,n
        v=gbs(:,j)*pi2a-v0
        v1=sqrt(dot_product(v,v))
        if(v1<1.d-8) then
          vc1(j)=amegaint/amega
        else
          do iatom=1,natom
            isort=is(iatom)
            phase=par*dot_product(v,tau(:,iatom))
            cc=-4.d0*pi*dcmplx(cos(phase),sin(phase))/amega
            a=v1*smt(isort)
            vc1(j)=vc1(j)+cc*(sin(a)-a*cos(a))/v1**3
          enddo
        endif
      enddo
      vc2=(0.d0,0.d0)
      i=indpw_gw(n2,iq)
      v0=gbs(:,i)*pi2a
      do j=1,n
        v=v0-gbs(:,j)*pi2a
        v1=sqrt(dot_product(v,v))
        if(v1<1.d-8) then
          vc2(j)=amegaint/amega
        else
          do iatom=1,natom
            isort=is(iatom)
            phase=par*dot_product(v,tau(:,iatom))
            cc=-4.d0*pi*dcmplx(cos(phase),sin(phase))/amega
            a=v1*smt(isort)
            vc2(j)=vc2(j)+cc*(sin(a)-a*cos(a))/v1**3
          enddo
        endif
      enddo
      vii=(0.d0,0.d0)
      j=1
      if(iq==1) j=2
      do i=j,n
        v=(pnt(:,iq)+gbs(:,i))*pi2a
        v1=dot_product(v,v)
        vii=vii+vc1(i)*vc2(i)/v1
      enddo
      vii=vii*8.d0*pi
      deallocate(vc1,vc2)
      end
