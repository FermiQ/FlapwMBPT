      subroutine k_for_band_plot
      use manager_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
      implicit none
	  integer :: i,i1,i2,k,j,i3,n,m,n1,n2,n3,n0
	  character*1 :: s(19),s1(19)
      real*8 :: r(3),d,v(3),v1(3),rv,sv,sv0,eta,ksi,ac,bc,dlt,mu,a,b,c,
     &          csa,sna,nu,cskg,eps,rho,cska,cskb,cskc,ab,ca2,cb,cb2,
     &          comp,ga,gb,gc,omega,phi,psi
	  integer, allocatable :: ind_symb(:)
	  real*8, allocatable :: v_symb(:,:)
	  allocate(v_symb(3,19),ind_symb(20))
	  eps=1.d-10
	  if(istruc==1) then   !!! simple cubic
	    n_symb=4
	    s(1:4)=(/'G','X','M','R'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	    v_symb(:,2)=(/0.5d0,0.d0,0.d0/)
	    v_symb(:,3)=(/0.5d0,0.5d0,0.d0/)
	    v_symb(:,4)=(/0.5d0,0.5d0,0.5d0/)
	    num_k_path=7
	    ind_symb(1:8)=(/1,2,3,1,4,2,4,3/)
	  else if(istruc==2) then   !!! Body-centered cubic
	    n_symb=4
	    s(1:4)=(/'G','H','P','N'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	    v_symb(:,2)=(/0.5d0,-0.5d0,0.5d0/)
	    v_symb(:,3)=(/0.25d0,0.25d0,0.25d0/)
	    v_symb(:,4)=(/0.d0,0.d0,0.5d0/)
	    num_k_path=7
	    ind_symb(1:8)=(/1,2,4,1,3,2,4,3/)
	  else if(istruc==3) then   !!! Face-centered cubic
	    n_symb=6
	    s(1:6)=(/'G','X','W','L','K','U'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	    v_symb(:,2)=(/0.5d0,0.d0,0.5d0/)
	    v_symb(:,3)=(/0.5d0,0.25d0,0.75d0/)
	    v_symb(:,4)=(/0.5d0,0.5d0,0.5d0/)
	    v_symb(:,5)=(/0.375d0,0.375d0,0.75d0/)
	    v_symb(:,6)=(/0.625d0,0.25d0,0.625d0/)
	    num_k_path=13
	    ind_symb(1:14)=(/1,5,3,2,1,4,6,1,3,6,2,5,4,3/)
	  else if(istruc==4) then   !!! Hexagonal
	    n_symb=6
	    s(1:6)=(/'G','K','M','A','H','L'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
        v_symb(:,2)=(/0.3333333333333333d0,0.33333333333333333d0,0.d0/)
	    v_symb(:,3)=(/0.d0,0.5d0,0.d0/)
	    v_symb(:,4)=(/0.d0,0.d0,0.5d0/)
        v_symb(:,5)=(/0.3333333333333333d0,0.33333333333333333d0,0.5d0/)
	    v_symb(:,6)=(/0.d0,0.5d0,0.5d0/)
	    num_k_path=12
	    ind_symb(1:13)=(/1,3,2,1,4,6,5,4,2,5,3,6,1/)
	  else if(istruc==5) then   !!! Simple tetragonal
	    n_symb=6
	    s(1:6)=(/'G','Z','A','M','R','X'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	    v_symb(:,2)=(/0.d0,0.d0,0.5d0/)
	    v_symb(:,3)=(/0.5d0,0.5d0,0.5d0/)
	    v_symb(:,4)=(/0.5d0,0.5d0,0.d0/)
	    v_symb(:,5)=(/0.5d0,0.d0,0.5d0/)
	    v_symb(:,6)=(/0.5d0,0.d0,0.d0/)
	    num_k_path=12
	    ind_symb(1:13)=(/1,6,4,1,2,5,3,2,6,5,4,3,1/)
	  else if(istruc==6) then   !!! Body-centered tetragonal  C < A
	    n_symb=7
	    s(1:7)=(/'G','M','N','P','X','Z','z'/)
	    eta=(1.d0+c_a**2)/4.d0
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	    v_symb(:,2)=(/-0.5d0,0.5d0,0.5d0/)
	    v_symb(:,3)=(/0.d0,0.5d0,0.d0/)
	    v_symb(:,4)=(/0.25d0,0.25d0,0.25d0/)
	    v_symb(:,5)=(/0.d0,0.d0,0.5d0/)
	    v_symb(1:2,6)=eta
        v_symb(3,6)=-eta
        v_symb(1,7)=-eta
        v_symb(2,7)=1.d0-eta
        v_symb(3,7)=eta
	    num_k_path=13
	    ind_symb(1:14)=(/1,5,2,1,6,4,3,7,2,4,5,7,4,1/)
      else if(istruc==7) then   !!! Body-centered tetragonal  C > A
        n_symb=9
        s(1:9)=(/'G','N','P','S','s','X','Y','y','Z'/)
        eta=(1.d0+c_a**2)/4.d0
        ksi=0.5d0/c_a**2
        v_symb(:,1)=(/0.d0,0.d0,0.d0/)
        v_symb(:,2)=(/0.d0,0.5d0,0.d0/)
        v_symb(:,3)=(/0.25d0,0.25d0,0.25d0/)
        v_symb(1,4)=-eta
        v_symb(2,4)=eta
        v_symb(3,4)=eta
        v_symb(1,5)=eta
        v_symb(2,5)=1.d0-eta
        v_symb(3,5)=-eta
        v_symb(:,6)=(/0.d0,0.d0,0.5d0/)
        v_symb(1,7)=-ksi
        v_symb(2,7)=ksi
        v_symb(3,7)=0.5d0
        v_symb(1:2,8)=-.5d0
        v_symb(3,8)=-ksi
        v_symb(:,9)=(/0.5d0,0.5d0,-0.5d0/)
        num_k_path=17
        ind_symb(1:18)=(/1,6,7,4,1,9,5,2,3,8,9,6,3,7,2,8,5,1/)
	  else if(istruc==8) then   !!! Simple orthorhombic
	    n_symb=8
	    s(1:8)=(/'G','R','S','T','U','X','Y','Z'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	    v_symb(:,2)=(/0.5d0,0.5d0,0.5d0/)
	    v_symb(:,3)=(/0.5d0,0.5d0,0.d0/)
	    v_symb(:,4)=(/0.d0,0.5d0,0.5d0/)
	    v_symb(:,5)=(/0.5d0,0.d0,0.5d0/)
	    v_symb(:,6)=(/0.5d0,0.d0,0.d0/)
	    v_symb(:,7)=(/0.d0,0.5d0,0.d0/)
	    v_symb(:,8)=(/0.d0,0.d0,0.5d0/)
	    num_k_path=16
	    ind_symb(1:17)=(/1,6,3,7,1,8,5,2,4,8,6,5,3,2,7,4,1/)
	  else if(istruc==9) then   !!! Body-centered orthorhombic
	    ac=1.d0/c_a
	    bc=b_a/c_a
        eta=(1.d0+bc**2)/4.d0
        ksi=(1.d0+ac**2)/4.d0
        dlt=(bc**2-ac**2)/4.d0
        mu=(ac**2+bc**2)/4.d0
	    n_symb=12
        s(1:12)=(/'G','L','l','R','S','T','W','X','x','Y','y','Z'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	    v_symb(1,2)=-mu
        v_symb(2,2)=mu
        v_symb(3,2)=0.5d0-dlt
        v_symb(1,3)=mu
        v_symb(2,3)=-mu
        v_symb(3,3)=0.5d0+dlt
	    v_symb(:,4)=(/0.d0,0.5d0,0.d0/)
	    v_symb(:,5)=(/0.5d0,0.d0,0.d0/)
	    v_symb(:,6)=(/0.d0,0.d0,0.5d0/)
	    v_symb(:,7)=(/0.25d0,0.25d0,0.25d0/)
        v_symb(1,8)=-ksi
        v_symb(2:3,8)=ksi
        v_symb(1,9)=ksi
        v_symb(2,9)=1.d0-ksi
        v_symb(3,9)=-ksi
        v_symb(1,10)=eta
        v_symb(2,10)=-eta
        v_symb(3,10)=eta
        v_symb(1,11)=1.d0-eta
        v_symb(2,11)=eta
        v_symb(3,11)=-eta
        v_symb(:,12)=(/0.5d0,0.5d0,-0.5d0/)
	    num_k_path=15
        ind_symb(1:16)=(/1,8,2,6,7,4,9,12,1,10,5,7,3,10,11,12/)
	  else if(istruc==10) then   !!! Face-centered orthorhombic
	    ca2=c_a**2
	    cb2=(c_a/b_a)**2+1.d0
	    if(ca2-cb2>eps) then
          ac=1.d0/c_a
          ab=1.d0/b_a
          eta=(1.d0+ab**2+ac**2)/4.d0
          ksi=(1.d0+ab**2-ac**2)/4.d0
	      n_symb=9
	      s(1:9)=(/'G','A','a','L','T','X','x','Y','Z'/)
	      v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(1,2)=0.5d0
          v_symb(2,2)=0.5d0+ksi
          v_symb(3,2)=ksi
          v_symb(1,3)=0.5d0
          v_symb(2,3)=0.5d0-ksi
          v_symb(3,3)=1.d0-ksi
	      v_symb(:,4)=(/0.5d0,0.5d0,0.5d0/)
	      v_symb(:,5)=(/1.d0,0.5d0,0.5d0/)
          v_symb(1,6)=0.d0
          v_symb(2:3,6)=eta
          v_symb(1,7)=1.d0
          v_symb(2:3,7)=1.d0-eta
	      v_symb(:,8)=(/0.5d0,0.d0,0.5d0/)
	      v_symb(:,9)=(/0.5d0,0.5d0,0.d0/)
	      num_k_path=14
	      ind_symb(1:15)=(/1,8,5,9,1,6,3,8,7,5,1,4,6,2,9/)
        else if(ca2-cb2<-eps) then
          ac=1.d0/c_a
          ab=1.d0/b_a
          cb=c_a/b_a
          eta=(1.d0+ab**2-ac**2)/4.d0
          phi=(1.d0+cb**2-c_a**2)/4.d0
          dlt=(1.d0+b_a**2-1.d0/cb**2)/4.d0
          n_symb=11
          s(1:11)=(/'G','C','c','D','d','L','H','h','X','Y','Z'/)
          v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(1,2)=0.5d0
          v_symb(2,2)=0.5d0-eta
          v_symb(3,2)=1.d0-eta
          v_symb(1,3)=0.5d0
          v_symb(2,3)=0.5d0+eta
          v_symb(3,3)=eta
          v_symb(1,4)=0.5d0-dlt
          v_symb(2,4)=0.5d0
          v_symb(3,4)=1.d0-dlt
          v_symb(1,5)=0.5d0+dlt
          v_symb(2,5)=0.5d0
          v_symb(3,5)=dlt
          v_symb(:,6)=(/0.5d0,0.5d0,0.5d0/)
          v_symb(1,7)=1.d0-phi
          v_symb(2,7)=0.5d0-phi
          v_symb(3,7)=0.5d0
          v_symb(1,8)=phi
          v_symb(2,8)=0.5d0+phi
          v_symb(3,8)=0.5d0
          v_symb(:,9)=(/0.d0,0.5d0,0.5d0/)
          v_symb(:,10)=(/0.5d0,0.d0,0.5d0/)
          v_symb(:,11)=(/0.5d0,0.5d0,0.d0/)
          num_k_path=17
          ind_symb(1:18)=(/1,10,2,4,9,1,11,5,7,2,1,6,9,8,3,11,10,7/)
        else if(abs(ca2-cb2)<eps) then
          ac=1.d0/c_a
          ab=1.d0/b_a
          eta=(1.d0+ab**2+ac**2)/4.d0
          ksi=(1.d0+ab**2-ac**2)/4.d0
          n_symb=9
          s(1:9)=(/'G','A','a','L','T','X','x','Y','Z'/)
          v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(1,2)=0.5d0
          v_symb(2,2)=0.5d0+ksi
          v_symb(3,2)=ksi
          v_symb(1,3)=0.5d0
          v_symb(2,3)=0.5d0-ksi
          v_symb(3,3)=1.d0-ksi
          v_symb(:,4)=(/0.5d0,0.5d0,0.5d0/)
          v_symb(:,5)=(/1.d0,0.5d0,0.5d0/)
          v_symb(1,6)=0.d0
          v_symb(2:3,6)=eta
          v_symb(1,7)=1.d0
          v_symb(2:3,7)=1.d0-eta
          v_symb(:,8)=(/0.5d0,0.d0,0.5d0/)
          v_symb(:,9)=(/0.5d0,0.5d0,0.d0/)
          num_k_path=15
          ind_symb(1:16)=(/1,8,5,9,1,6,3,8,6,2,9,3,5,2,1,4/)
        endif
      else if(istruc==11) then   !!! Base-centered orthorhombic
        ab=1.d0/b_a
        ksi=(1.d0+ab**2)/4.d0
        n_symb=10
        s(1:10)=(/'G','A','a','R','S','T','X','x','Y','Z'/)
        v_symb(:,1)=(/0.d0,0.d0,0.d0/)
        v_symb(1:2,2)=ksi
        v_symb(3,2)=0.5d0
        v_symb(1,3)=-ksi
        v_symb(2,3)=1.d0-ksi
        v_symb(3,3)=0.5d0
        v_symb(:,4)=(/0.d0,0.5d0,0.5d0/)
        v_symb(:,5)=(/0.d0,0.5d0,0.d0/)
        v_symb(:,6)=(/-0.5d0,0.5d0,0.5d0/)
        v_symb(1:2,7)=ksi
        v_symb(3,7)=0.d0
        v_symb(1,8)=-ksi
        v_symb(2,8)=1.d0-ksi
        v_symb(3,8)=0.d0
        v_symb(:,9)=(/-0.5d0,0.5d0,0.d0/)
        v_symb(:,10)=(/0.d0,0.d0,0.5d0/)
        num_k_path=17
        ind_symb(1:18)=(/1,7,5,4,2,10,1,9,8,3,6,9,10,6,8,5,1,4/)
	  else if(istruc==12) then   !!! Simple monoclinic
c     C axis is unique along Z; A axis is along X
c     Lattice:  A=(a,0,0)
c               B=(b*cos(gamma),b*sin(gamma),0)
c               C=(0,0,c)
        a=rbas(1,1)
        b=sqrt(rbas(1,2)**2+rbas(2,2)**2)
        c=rbas(3,3)
        csa=rbas(1,2)/b
        sna=rbas(2,2)/b
        eta=(1.d0+a*csa/b)/2.d0/sna**2
        nu=0.5d0+eta*b*csa/a
	    n_symb=16
        s(1:16)=(/'G','A','C','D','d','E','H','h','h','M','m','m','X',
     &            'Y','y','Z'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	    v_symb(:,2)=(/0.5d0,0.0d0,0.5d0/)
	    v_symb(:,3)=(/0.5d0,0.5d0,0.d0/)
	    v_symb(:,4)=(/0.d0,0.5d0,0.5d0/)
	    v_symb(:,5)=(/0.d0,-0.5d0,0.5d0/)
	    v_symb(:,6)=(/0.5d0,0.5d0,0.5d0/)
	    v_symb(1,7)=eta
        v_symb(2,7)=1.d0-nu
        v_symb(3,7)=0.d0
        v_symb(1,8)=1.d0-eta
        v_symb(2,8)=nu
        v_symb(3,8)=0.d0
        v_symb(1,9)=eta
        v_symb(2,9)=-nu
        v_symb(3,9)=0.d0
        v_symb(1,10)=eta
        v_symb(2,10)=1.d0-nu
        v_symb(3,10)=0.5d0
        v_symb(1,11)=1.d0-eta
        v_symb(2,11)=nu
        v_symb(3,11)=0.5d0
        v_symb(1,12)=eta
        v_symb(2,12)=-nu
        v_symb(3,12)=0.5d0
	    v_symb(:,13)=(/0.5d0,0.d0,0.d0/)
        v_symb(:,14)=(/0.d0,0.5d0,0.d0/)
        v_symb(:,15)=(/0.d0,-0.5d0,0.d0/)
        v_symb(:,16)=(/0.d0,0.d0,0.5d0/)
        num_k_path=15
        ind_symb(1:16)=(/1,14,7,3,6,11,2,13,8,1,3,10,4,16,14,4/)
      else if(istruc==121) then   !!! Simple monoclinic
c     a,b<c;alpha<90;beta,gamma=90
c     Lattice:  A=(a,0,0)
c               B=(0,b,0)
c               C=(0,c*cos(alpha),c*sin(alpha))
        a=rbas(1,1)
        b=rbas(2,2)
        c=sqrt(rbas(2,3)**2+rbas(3,3)**2)
        csa=rbas(2,3)/c
        sna=rbas(3,3)/c
        eta=(1.d0-b*csa/c)/2.d0/sna**2
        nu=0.5d0-eta*c*csa/b
        n_symb=18
        s(1:18)=(/'G','A','C','D','d','E','H','h','h','M','m','m','X',
     &            'Y','y','Z','Q','R'/)
        v_symb(:,1)=(/0.d0,0.d0,0.d0/)
        v_symb(:,2)=(/0.5d0,0.5d0,0.d0/)
        v_symb(:,3)=(/0.d0,0.5d0,0.5d0/)
        v_symb(:,4)=(/0.5d0,0.d0,0.d0/)
        v_symb(:,5)=(/0.5d0,0.d0,-0.5d0/)
        v_symb(:,6)=(/0.5d0,0.5d0,0.5d0/)
        v_symb(1,7)=0.d0
        v_symb(2,7)=eta
        v_symb(3,7)=1.d0-nu
        v_symb(1,8)=0.d0
        v_symb(2,8)=1.d0-eta
        v_symb(3,8)=nu
        v_symb(1,9)=0.d0
        v_symb(2,9)=eta
        v_symb(3,9)=-nu
        v_symb(1,10)=0.5d0
        v_symb(2,10)=eta
        v_symb(3,10)=1.d0-nu
        v_symb(1,11)=0.5d0
        v_symb(2,11)=1.d0-eta
        v_symb(3,11)=nu
        v_symb(1,12)=0.5d0
        v_symb(2,12)=eta
        v_symb(3,12)=-nu
        v_symb(:,13)=(/0.d0,0.5d0,0.d0/)
        v_symb(:,14)=(/0.d0,0.d0,0.5d0/)
        v_symb(:,15)=(/0.d0,0.d0,-0.5d0/)
        v_symb(:,16)=(/0.5d0,0.d0,0.d0/)
        v_symb(:,17)=(/0.d0,-0.25d0,0.25d0/)
        v_symb(:,18)=(/0.5d0,-0.25d0,0.25d0/)
        num_k_path=13
        ind_symb(1:14)=(/1,14,7,3,6,11,2,13,8,1,17,14,4,16/)
c        num_k_path=15
c        ind_symb(1:16)=(/1,14,7,3,6,11,2,13,8,1,3,10,4,16,14,4/)
c ----------------- Original ---------------------------------------
c       num_k_path=15
c       ind_symb(1:16)=(/1,14,7,3,6,11,2,13,8,1,10,4,16,14,4,1/)
c ----------------- Paper + Q, R added -----------------------------
c        num_k_path=17
c        ind_symb(1:18)=(/1,14,7,3,6,11,2,13,8,10,4,16,18,17,4,14,17,1/)
      else if(istruc==13) then   !!! Base-centered monoclinic
c     a,b<c;alpha<90;beta,gamma=90
c     Lattice:  A=(a/2,b/2,0)
c               B=(-a/2,b/2,0)
c               C=(0,c*cos(alpha),c*sin(alpha))
        a=2.d0*rbas(1,1)
        b=2.d0*rbas(2,1)
        c=sqrt(rbas(2,3)**2+rbas(3,3)**2)
        csa=rbas(2,3)/c
        sna=rbas(3,3)/c
        ga=sqrt(dot_product(gbas(:,1),gbas(:,1)))
        gb=sqrt(dot_product(gbas(:,2),gbas(:,2)))
        cskg=dot_product(gbas(:,1),gbas(:,2))/ga/gb
        comp=b*csa/c+(b*sna/a)**2
        if(cskg<-eps) then
          ksi=(2.d0-b*csa/c)/4.d0/sna**2
          eta=0.5d0+2.d0*ksi*c*csa/b
          psi=0.75d0-a*a/4.d0/b/b/sna**2
          phi=psi+(0.75d0-psi)*b*csa/c
          n_symb=17
          s(1:17)=(/'G','N','n','F','f','s','t','I','i','L','M','X','x',
     &              'j','Y','y','Z'/)
          v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(:,2)=(/0.5d0,0.d0,0.d0/)
          v_symb(:,3)=(/0.d0,-0.5d0,0.d0/)
          v_symb(1:2,4)=1.d0-ksi
          v_symb(3,4)=1.d0-eta
          v_symb(1:2,5)=ksi
          v_symb(3,5)=eta
          v_symb(1:2,6)=-ksi
          v_symb(3,6)=1.d0-eta
          v_symb(1,7)=1.d0-ksi
          v_symb(2,7)=-ksi
          v_symb(3,7)=1.d0-eta
          v_symb(1,8)=phi
          v_symb(2,8)=1.d0-phi
          v_symb(3,8)=0.5d0
          v_symb(1,9)=1.d0-phi
          v_symb(2,9)=phi-1.d0
          v_symb(3,9)=0.5d0
          v_symb(:,10)=(/0.5d0,0.5d0,0.5d0/)
          v_symb(:,11)=(/0.5d0,0.d0,0.5d0/)
          v_symb(1,12)=1.d0-psi
          v_symb(2,12)=psi-1.d0
          v_symb(3,12)=0.d0
          v_symb(1,13)=psi
          v_symb(2,13)=1.d0-psi
          v_symb(3,13)=0.d0
          v_symb(1,14)=psi-1.d0
          v_symb(2,14)=-psi
          v_symb(3,14)=0.d0
          v_symb(:,15)=(/0.5d0,0.5d0,0.d0/)
          v_symb(:,16)=(/-0.5d0,0.5d0,0.d0/)
          v_symb(:,17)=(/0.d0,0.d0,0.5d0/)
          num_k_path=14
          ind_symb(1:15)=(/1,15,4,10,8,9,17,5,15,13,2,1,12,11,1/)
        else if(abs(cskg)<eps) then
          ksi=(2.d0-b*csa/c)/4.d0/sna**2
          eta=0.5d0+2.d0*ksi*c*csa/b
          psi=0.75d0-a*a/4.d0/b/b/sna**2
          phi=psi+(0.75d0-psi)*b*csa/c
          n_symb=17
          s(1:17)=(/'G','N','n','F','f','s','t','I','i','L','M','X','x',
     &              'j','Y','y','Z'/)
          v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(:,2)=(/0.5d0,0.d0,0.d0/)
          v_symb(:,3)=(/0.d0,-0.5d0,0.d0/)
          v_symb(1:2,4)=1.d0-ksi
          v_symb(3,4)=1.d0-eta
          v_symb(1:2,5)=ksi
          v_symb(3,5)=eta
          v_symb(1:2,6)=-ksi
          v_symb(3,6)=1.d0-eta
          v_symb(1,7)=1.d0-ksi
          v_symb(2,7)=-ksi
          v_symb(3,7)=1.d0-eta
          v_symb(1,8)=phi
          v_symb(2,8)=1.d0-phi
          v_symb(3,8)=0.5d0
          v_symb(1,9)=1.d0-phi
          v_symb(2,9)=phi-1.d0
          v_symb(3,9)=0.5d0
          v_symb(:,10)=(/0.5d0,0.5d0,0.5d0/)
          v_symb(:,11)=(/0.5d0,0.d0,0.5d0/)
          v_symb(1,12)=1.d0-psi
          v_symb(2,12)=psi-1.d0
          v_symb(3,12)=0.d0
          v_symb(1,13)=psi
          v_symb(2,13)=1.d0-psi
          v_symb(3,13)=0.d0
          v_symb(1,14)=psi-1.d0
          v_symb(2,14)=-psi
          v_symb(3,14)=0.d0
          v_symb(:,15)=(/0.5d0,0.5d0,0.d0/)
          v_symb(:,16)=(/-0.5d0,0.5d0,0.d0/)
          v_symb(:,17)=(/0.d0,0.d0,0.5d0/)
          num_k_path=10
          ind_symb(1:11)=(/1,15,4,10,8,9,17,5,2,1,11/)
        else if(cskg>eps) then
          if(comp-1.d0<-eps) then
            mu=(1.d0+(b/a)**2)/4.d0
            dlt=b*c*csa/2.d0/a/a
            ksi=mu-0.25d0+(1.d0-b*csa/c)/4.d0/sna**2
            eta=0.5d0+2.d0*ksi*c*csa/b
            phi=1.d0+ksi-2.d0*mu
            psi=eta-2.d0*dlt
            n_symb=17
            s(1:17)=(/'G','F','f','i','H','h','j','I','M','N','n','X',
     &                'Y','y','k','l','Z'/)
            v_symb(:,1)=(/0.d0,0.d0,0.d0/)
            v_symb(1:2,2)=1.d0-phi
            v_symb(3,2)=1.d0-psi
            v_symb(1,3)=phi
            v_symb(2,3)=phi-1.d0
            v_symb(3,3)=psi
            v_symb(1,4)=1.d0-phi
            v_symb(2,4)=-phi
            v_symb(3,4)=1.d0-psi
            v_symb(1:2,5)=ksi
            v_symb(3,5)=eta
            v_symb(1,6)=1.d0-ksi
            v_symb(2,6)=-ksi
            v_symb(3,6)=1.d0-eta
            v_symb(1:2,7)=-ksi
            v_symb(3,7)=1.d0-eta
            v_symb(:,8)=(/0.5d0,-0.5d0,0.5d0/)
            v_symb(:,9)=(/0.5d0,0.d0,0.5d0/)
            v_symb(:,10)=(/0.5d0,0.d0,0.d0/)
            v_symb(:,11)=(/0.d0,-0.5d0,0.d0/)
            v_symb(:,12)=(/0.5d0,-0.5d0,0.d0/)
            v_symb(1:2,13)=mu
            v_symb(3,13)=dlt
            v_symb(1,14)=1.d0-mu
            v_symb(2,14)=-mu
            v_symb(3,14)=-dlt
            v_symb(1:2,15)=-mu
            v_symb(3,15)=-dlt
            v_symb(1,16)=mu
            v_symb(2,16)=mu-1.d0
            v_symb(3,16)=dlt
            v_symb(:,17)=(/0.d0,0.d0,0.5d0/)
            num_k_path=13
            ind_symb(1:14)=(/1,13,2,5,17,8,3,6,14,12,1,10,9,1/)
          else if(abs(comp-1.d0)<eps) then
            mu=(1.d0+(b/a)**2)/4.d0
            dlt=b*c*csa/2.d0/a/a
            ksi=mu-0.25d0+(1.d0-b*csa/c)/4.d0/sna**2
            eta=0.5d0+2.d0*ksi*c*csa/b
            phi=1.d0+ksi-2.d0*mu
            psi=eta-2.d0*dlt
            n_symb=17
            s(1:17)=(/'G','F','f','i','H','h','j','I','M','N','n','X',
     &                'Y','y','k','l','Z'/)
            v_symb(:,1)=(/0.d0,0.d0,0.d0/)
            v_symb(1:2,2)=1.d0-phi
            v_symb(3,2)=1.d0-psi
            v_symb(1,3)=phi
            v_symb(2,3)=phi-1.d0
            v_symb(3,3)=psi
            v_symb(1,4)=1.d0-phi
            v_symb(2,4)=-phi
            v_symb(3,4)=1.d0-psi
            v_symb(1:2,5)=ksi
            v_symb(3,5)=eta
            v_symb(1,6)=1.d0-ksi
            v_symb(2,6)=-ksi
            v_symb(3,6)=1.d0-eta
            v_symb(1:2,7)=-ksi
            v_symb(3,7)=1.d0-eta
            v_symb(:,8)=(/0.5d0,-0.5d0,0.5d0/)
            v_symb(:,9)=(/0.5d0,0.d0,0.5d0/)
            v_symb(:,10)=(/0.5d0,0.d0,0.d0/)
            v_symb(:,11)=(/0.d0,-0.5d0,0.d0/)
            v_symb(:,12)=(/0.5d0,-0.5d0,0.d0/)
            v_symb(1:2,13)=mu
            v_symb(3,13)=dlt
            v_symb(1,14)=1.d0-mu
            v_symb(2,14)=-mu
            v_symb(3,14)=-dlt
            v_symb(1:2,15)=-mu
            v_symb(3,15)=-dlt
            v_symb(1,16)=mu
            v_symb(2,16)=mu-1.d0
            v_symb(3,16)=dlt
            v_symb(:,17)=(/0.d0,0.d0,0.5d0/)
            num_k_path=12
            ind_symb(1:13)=(/1,13,2,5,17,8,6,14,12,1,10,9,1/)
          else if(comp-1.d0>eps) then
            ksi=((b/a)**2+(1.d0-b*csa/c)/sna**2)/4.d0
            eta=0.5d0+2.d0*ksi*c*csa/b
            mu=eta/2.d0+(b/2.d0/a)**2-b*c*csa/2.d0/a/a
            nu=2.d0*mu-ksi
            omega=(4.d0*nu-1.d0-(b*sna/a)**2)*c/2.d0/b/csa
            dlt=ksi*c*csa/b+omega/2.d0-0.25d0
            rho=1.d0-ksi*(a/b)**2
            n_symb=18
            s(1:18)=(/'G','F','f','i','H','h','j','I','k','M','N','n',
     &                'X','Y','y','l','m','Z'/)
            v_symb(:,1)=(/0.d0,0.d0,0.d0/)
            v_symb(1:2,2)=nu
            v_symb(3,2)=omega
            v_symb(1:2,3)=1.d0-nu
            v_symb(3,3)=1.d0-omega
            v_symb(1,4)=nu
            v_symb(2,4)=nu-1.d0
            v_symb(3,4)=omega
            v_symb(1:2,5)=ksi
            v_symb(3,5)=eta
            v_symb(1,6)=1.d0-ksi
            v_symb(2,6)=-ksi
            v_symb(3,6)=1.d0-eta
            v_symb(1:2,7)=-ksi
            v_symb(3,7)=1.d0-eta
            v_symb(1,8)=rho
            v_symb(2,8)=1.d0-rho
            v_symb(3,8)=0.5d0
            v_symb(1,9)=1.d0-rho
            v_symb(2,9)=rho-1.d0
            v_symb(3,9)=0.5d0
            v_symb(:,10)=(/0.5d0,0.5d0,0.5d0/)
            v_symb(:,11)=(/0.5d0,0.d0,0.5d0/)
            v_symb(:,12)=(/0.5d0,0.d0,0.d0/)
            v_symb(:,13)=(/0.d0,-0.5d0,0.d0/)
            v_symb(:,14)=(/0.5d0,-0.5d0,0.d0/)
            v_symb(1:2,15)=mu
            v_symb(3,15)=dlt
            v_symb(1,16)=1.d0-mu
            v_symb(2,16)=-mu
            v_symb(3,16)=-dlt
            v_symb(1:2,17)=-mu
            v_symb(3,17)=-dlt
            v_symb(1,18)=mu
            v_symb(2,18)=mu-1.d0
            v_symb(3,18)=dlt
            v_symb(:,19)=(/0.d0,0.d0,0.5d0/)
            num_k_path=15
            ind_symb(1:16)=(/1,15,2,10,8,9,19,5,3,6,16,14,1,12,11,1/)
          endif
        endif
	  else if(istruc==14) then   !!! Rombohedral
c     Lattice:
c     A=(a*cos(alpha/2),-a*sin(alpha/2),0)
c     B=(a*cos(alpha/2),-a*sin(alpha/2),0)
c     C=.....
        a=sqrt(dot_product(rbas(:,1),rbas(:,1)))
        csa=dot_product(rbas(:,1),rbas(:,2))/a/a
        if(csa>0.d0) then
          eta=(1.d0+4.d0*csa)/(2.d0+4.d0*csa)
          nu=0.75d0-eta/2.d0
	      n_symb=12
	      s(1:12)=(/'G','B','b','F','L','l','P','p','i','Q','X','Z'/)
	      v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(1,2)=eta
          v_symb(2,2)=0.5d0
          v_symb(3,2)=1.d0-eta
          v_symb(1,3)=0.5d0
          v_symb(2,3)=1.d0-eta
          v_symb(3,3)=eta-1.d0
	      v_symb(:,4)=(/0.5d0,0.5d0,0.d0/)
          v_symb(:,5)=(/0.5d0,0.d0,0.d0/)
          v_symb(:,6)=(/0.d0,0.d0,-0.5d0/)
          v_symb(1,7)=eta
          v_symb(2:3,7)=nu
          v_symb(1:2,8)=1.d0-nu
          v_symb(3,8)=1.d0-eta
          v_symb(1:2,9)=nu
          v_symb(3,9)=eta-1.d0
          v_symb(1,10)=1.d0-nu
          v_symb(2:3,10)=nu
          v_symb(1,11)=nu
          v_symb(2,11)=0.d0
          v_symb(3,11)=-nu
          v_symb(:,12)=(/0.5d0,0.5d0,0.5d0/)
	      num_k_path=12
	      ind_symb(1:13)=(/1,5,3,2,12,1,11,10,4,8,12,5,7/)
        else if(csa<0.d0) then
          eta=(1.d0+csa)/2.d0/(1.d0-csa)
          nu=0.75d0-eta/2.d0
          n_symb=8
          s(1:8)=(/'G','F','L','P','p','Q','q','Z'/)
          v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(:,2)=(/0.5d0,-0.5d0,0.d0/)
          v_symb(:,3)=(/0.5d0,0.d0,0.d0/)
          v_symb(1,4)=1.d0-nu
          v_symb(2,4)=-nu
          v_symb(3,4)=1.d0-nu
          v_symb(1,5)=nu
          v_symb(2:3,5)=nu-1.d0
          v_symb(:,6)=eta
          v_symb(1,7)=1.d0-eta
          v_symb(2:3,7)=-eta
          v_symb(:,8)=(/0.5d0,-0.5d0,0.5d0/)
          num_k_path=9
          ind_symb(1:10)=(/1,4,8,6,1,2,5,7,3,8/)
        endif
      else if(istruc==15) then   !!! Triclinic
        a=2.d0*rbas(1,1)
        b=2.d0*rbas(2,1)
        c=sqrt(rbas(2,3)**2+rbas(3,3)**2)
        csa=rbas(2,3)/c
        sna=rbas(3,3)/c
        ga=sqrt(dot_product(gbas(:,1),gbas(:,1)))
        gb=sqrt(dot_product(gbas(:,2),gbas(:,2)))
        gc=sqrt(dot_product(gbas(:,3),gbas(:,3)))
        cska=dot_product(gbas(:,2),gbas(:,3))/gb/gc
        cskb=dot_product(gbas(:,1),gbas(:,3))/ga/gc
        cskc=dot_product(gbas(:,1),gbas(:,2))/ga/gb
        if(cska<-eps.and.cskb<-eps) then
          n_symb=8
          s(1:8)=(/'G','L','M','N','R','X','Y','Z'/)
          v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(:,2)=(/0.5d0,0.5d0,0.d0/)
          v_symb(:,3)=(/0.d0,0.5d0,0.5d0/)
          v_symb(:,4)=(/0.5d0,0.d0,0.5d0/)
          v_symb(:,5)=(/0.5d0,0.5d0,0.5d0/)
          v_symb(:,6)=(/0.5d0,0.d0,0.d0/)
          v_symb(:,7)=(/0.d0,0.5d0,0.d0/)
          v_symb(:,8)=(/0.d0,0.d0,0.5d0/)
          num_k_path=10
          ind_symb(1:11)=(/6,1,7,2,1,8,4,1,3,5,1/)
        else if(cska>eps.and.cskb>eps) then
          n_symb=8
          s(1:8)=(/'G','L','M','N','R','X','Y','Z'/)
          v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(:,2)=(/0.5d0,-0.5d0,0.d0/)
          v_symb(:,3)=(/0.d0,0.d0,0.5d0/)
          v_symb(:,4)=(/-0.5d0,-0.5d0,0.5d0/)
          v_symb(:,5)=(/0.d0,-0.5d0,0.5d0/)
          v_symb(:,6)=(/0.d0,-0.5d0,0.d0/)
          v_symb(:,7)=(/0.5d0,0.d0,0.d0/)
          v_symb(:,8)=(/-0.5d0,0.d0,0.5d0/)
          num_k_path=10
          ind_symb(1:11)=(/6,1,7,2,1,8,4,1,3,5,1/)
        else
          n_symb=8
          s(1:8)=(/'G','L','M','N','R','X','Y','Z'/)
          v_symb(:,1)=(/0.d0,0.d0,0.d0/)
          v_symb(:,2)=(/0.5d0,0.5d0,0.d0/)
          v_symb(:,3)=(/0.d0,0.5d0,0.5d0/)
          v_symb(:,4)=(/0.5d0,0.d0,0.5d0/)
          v_symb(:,5)=(/0.5d0,0.5d0,0.5d0/)
          v_symb(:,6)=(/0.5d0,0.d0,0.d0/)
          v_symb(:,7)=(/0.d0,0.5d0,0.d0/)
          v_symb(:,8)=(/0.d0,0.d0,0.5d0/)
          num_k_path=10
          ind_symb(1:11)=(/6,1,7,2,1,8,4,1,3,5,1/)
        endif
	  else if(istruc==16) then   !!! Hexagonal_two_dimensional
	    n_symb=3
	    s(1:3)=(/'G','K','M'/)
	    v_symb(:,1)=(/0.d0,0.d0,0.d0/)
        v_symb(:,2)=(/0.3333333333333333d0,0.33333333333333333d0,0.d0/)
	    v_symb(:,3)=(/0.d0,0.5d0,0.d0/)
	    num_k_path=4
	    ind_symb(1:5)=(/1,3,1,2,1/)
  	  endif
c -------------------------------------------------------------
	  do i=1,n_symb
	    r=v_symb(1,i)*gbas(:,1)+v_symb(2,i)*gbas(:,2)
     &	                       +v_symb(3,i)*gbas(:,3)
	    v_symb(:,i)=r
	  enddo
	  allocate(k_label(2,num_k_path))
	  allocate(kpath_len(num_k_path))
	  do i=1,num_k_path
	    i1=ind_symb(i)
	    i2=ind_symb(i+1)
	    k_label(1,i)=s(i1)
	    k_label(2,i)=s(i2)
	    r=v_symb(:,i2)-v_symb(:,i1)
	    d=sqrt(dot_product(r,r))
	    kpath_len(i)=d
	  enddo
c ---------------------------------------------------------------
	  num_k_plot=num_k_path*n_k_div+1
	  allocate(k_dist(num_k_plot))
	  k_dist(1)=0.d0
	  allocate(k_plot(3,num_k_plot))
	  k_plot(:,1)=v_symb(:,1)
	  k=1
	  do i=1,num_k_path
	    i1=ind_symb(i)
	    i2=ind_symb(i+1)
	    r=v_symb(:,i2)-v_symb(:,i1)
	    d=sqrt(dot_product(r,r))
	    r=r/n_k_div
	    d=d/n_k_div
	    do j=1,n_k_div
	      k=k+1
	      k_plot(:,k)=v_symb(:,i1)+j*r
	      k_dist(k)=k_dist(k-1)+d
	    enddo
	  enddo
	  allocate(ndim_kplot(nproc_k))
	  allocate(n_mpi_kplot(nproc_k))
	  call size_shift_par(num_k_plot,nproc_k,ndim_kplot,n_mpi_kplot)
c	 We define coordinates of points in 1 BZ ------------
	  allocate(pnt_1bz(3,product(ndiv+1)))
	  k=0
      do i3=0,ndiv(3)
        do i2=0,ndiv(2)
          do i1=0,ndiv(1)
            k=k+1
	        pnt_1bz(:,k)=i1*qb0(:,1)+i2*qb0(:,2)+i3*qb0(:,3)
	      enddo
	    enddo
	  enddo
c ------- We shift K_PLOT points into 0 BZ -----------------------
	  do k=1,num_k_plot
        call zone0_cart(k_plot(1,k),gbas,k_plot(1,k))
	  enddo
c ---- Associate symmetry label with k-points (if applicable) ----
      allocate(sym_label(npnt))
      sym_label=' '
      do i=1,n_symb
	    call zone1_number(v_symb(1,i),rb0,ndiv,i1)
	    if(i1>0) then
	      i1=index_k1(i1)
	      i1=i_kref(i1)
	      sym_label(i1)=s(i)
	    endif
	  enddo
c ------ Entry points for ksym in the NPNT - list -------------------
      allocate(ksym_in_npnt(n_symb))
      ksym_in_npnt=0
      allocate(k_npnt_in_sym(npnt))
      k_npnt_in_sym=0
      n=0
      m=0
      do i=1,n_symb
        i3=0
        do i2=1,n
          if(s(i)==s1(i2)) i3=1
        enddo
        if(i3==0) then
	      n=n+1
	      s1(n)=s(i)
	      call zone1_number(v_symb(1,i),rb0,ndiv,i1)
	      if(i1>0) then
	        i1=index_k1(i1)
	        m=m+1
	        ksym_in_npnt(m)=i_kref(i1)
	      endif
	    endif
	  enddo
	  n_sym=m
	  allocate(q_sym(3,n_sym))
      allocate(sym0_label(n_sym))
	  do k=1,n_sym
	    m=ksym_in_npnt(k)
	    q_sym(:,k)=pnt(:,m)
	    k_npnt_in_sym(m)=k
	    sym0_label(k)=sym_label(m)
	  enddo
	  if(maswrk) then
	    write(iun,*)' The following SYM K-points are for plotting'
	    do k=1,n_sym
          write(iun,'(3(1x,f8.4),a6,1x,f8.4)')(q_sym(i,k),i=1,3),
     &     '!K! = ',sqrt(dot_product(q_sym(:,k),q_sym(:,k)))
	    enddo
	  endif
c ------ Entry points for kline in the NPNT - list -------------------
      allocate(kline_in_npnt(max(n_symb,maxval(ndiv))))
      kline_in_npnt=0
      allocate(k_npnt_in_line(npnt))
      k_npnt_in_line=0
      if(k_line=='ALL') then
        m=npnt
        do i=1,npnt
          kline_in_npnt(i)=i
        enddo
      else 
        if(k_line=='100') then
          r=0.d0
          r(1)=1.d0
        else if(k_line=='010') then
          r=0.d0
          r(2)=1.d0
        else if(k_line=='001') then
          r=0.d0
          r(3)=1.d0
        else if(k_line=='110') then
          r=1.d0
          r(3)=0.d0
        else if(k_line=='111') then
          r=1.d0
        endif
        r=r/sqrt(dot_product(r,r))
        n1=ndiv(1)
        n2=ndiv(2)
        n3=ndiv(3)
        m=0
        do i3=-n3,n3
          do i2=-n2,n2
            do i1=-n1,n1
              v(1)=i1*qb0(1,1)+i2*qb0(1,2)+i3*qb0(1,3)
              v(2)=i1*qb0(2,1)+i2*qb0(2,2)+i3*qb0(2,3)
              v(3)=i1*qb0(3,1)+i2*qb0(3,2)+i3*qb0(3,3)
              if(dot_product(v,v)<1.d-10) then
                rv=0.d0
              else
                v1=v/sqrt(dot_product(v,v))
                rv=dot_product(r,v1)-1.d0
              endif
              if(abs(rv)<1.d-12) then
	            call zone1_number(v,rb0,ndiv,n)
	            n=i_kref(index_k1(n))
	            j=m+1
	            sv=dot_product(pnt(:,n),pnt(:,n))
	            do i=1,m
                  n0=kline_in_npnt(i)
	              sv0=dot_product(pnt(:,n0),pnt(:,n0))
	              if(abs(sv-sv0)<1.d-8) j=0
	            enddo
	            if(j/=0) then
	              m=m+1
	              kline_in_npnt(m)=n
	            endif
	          endif
	        enddo
	      enddo
	    enddo
	    do i1=1,m-1
	      i=kline_in_npnt(i1)
	      rv=dot_product(pnt(:,i),pnt(:,i))
	      i3=i1
	      do i2=i1+1,m
	        j=kline_in_npnt(i2)
	        d=dot_product(pnt(:,j),pnt(:,j))
	        if(d<rv) then
	          rv=d
	          i3=i2
	          n=j
	        endif
	      enddo
	      if(i3/=i1) then
	        kline_in_npnt(i1)=n
	        kline_in_npnt(i3)=i
	      endif
	    enddo
	  endif
	  n_line=m
	  allocate(q_line(3,n_line))
      if(l_double_freq) allocate(k_npnt_in_line_c(npnt_c))
      if(l_double_freq) then
        nc_line=0
        k_npnt_in_line_c=0
      endif
	  do k=1,n_line
	    q_line(:,k)=pnt(:,kline_in_npnt(k))
        i=kline_in_npnt(k)
	    k_npnt_in_line(i)=k
        if(l_double_freq) then
          j=k_c_from_a(i)
          if(j/=0) then
            nc_line=nc_line+1
            k_npnt_in_line_c(j)=nc_line
          endif
        endif
      enddo
      if(l_double_freq) then
        allocate(kline_in_npnt_c(nc_line))
	    do k=1,npnt_c
          i=k_npnt_in_line_c(k)
          if(i==0) cycle
          kline_in_npnt_c(i)=k
        enddo
      endif
	  if(maswrk) then
	    write(iun,*)' The following LINE K-points are for plotting'
	    do k=1,n_line
          write(iun,'(3(1x,f8.4),a6,1x,f8.4)')(q_line(i,k),i=1,3),
     &     '!K! = ',sqrt(dot_product(q_line(:,k),q_line(:,k)))
	    enddo
      endif
      end
