      subroutine p_head_num(ph,pha,ph1)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: ph(n_k_head_0,ndim3_nu)
      real*8, intent(out) :: pha(6,ndim3_nu),ph1(n_phead_dir,ndim3_nu)
      integer :: i,ind_nu,info,iq,k,k1
      real*8 :: pi2a,v(3),q2,b,c
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: a(:,:),ve(:,:),ph0(:,:)
      pi2a=(pi+pi)/par
c ------ Evaluate P_head for each direction ------------------------
      allocate(ph0(n_phead_dir,ndim3_nu))
      if(n_phead_pnt==1) then
        do i=1,n_phead_dir
          k=khead(1,i)
          iq=khead0(k)
          q2=pi2a**2*dot_product(pnt(:,iq),pnt(:,iq))
          do ind_nu=1,ndim3_nu
            ph0(i,ind_nu)=(ph(k,ind_nu)-ph(1,ind_nu))/q2
          enddo
        enddo
      else  if(n_phead_pnt==2) then
        do i=1,n_phead_dir
          k=khead(1,i)
          k1=khead(2,i)
          iq=khead0(k)
          q2=pi2a**2*dot_product(pnt(:,iq),pnt(:,iq))
          do ind_nu=1,ndim3_nu
            b=ph(k,ind_nu)-ph(1,ind_nu)
            c=ph(k1,ind_nu)-ph(1,ind_nu)
            ph0(i,ind_nu)=(8*b-c)/4.d0/q2
          enddo
        enddo
      endif
      ph1=ph0
c ------------------------------------------------------------------
      pha=0.d0
      if(istruc<=3) then
        do ind_nu=1,ndim3_nu
          i=ind_tens33(1,1)
          pha(i,ind_nu)=ph0(1,ind_nu)
          i=ind_tens33(2,2)
          pha(i,ind_nu)=ph0(1,ind_nu)
          i=ind_tens33(3,3)
          pha(i,ind_nu)=ph0(1,ind_nu)
        enddo
      else if((istruc>3.and.istruc<8).or.istruc==14) then
        do ind_nu=1,ndim3_nu
          i=ind_tens33(1,1)
          pha(i,ind_nu)=ph0(1,ind_nu)
          i=ind_tens33(2,2)
          pha(i,ind_nu)=ph0(1,ind_nu)
          i=ind_tens33(3,3)
          pha(i,ind_nu)=ph0(2,ind_nu)
        enddo
      else if(istruc>7.and.istruc<12) then
        do ind_nu=1,ndim3_nu
          i=ind_tens33(1,1)
          pha(i,ind_nu)=ph0(1,ind_nu)
          i=ind_tens33(2,2)
          pha(i,ind_nu)=ph0(2,ind_nu)
          i=ind_tens33(3,3)
          pha(i,ind_nu)=ph0(3,ind_nu)
        enddo
      else if(istruc==12) then
        allocate(ve(3,3))
        do i=1,3
          if(i==1) v=qb0(:,1)
          if(i==2) v=qb0(:,1)-qb0(:,2)
          if(i==3) v=qb0(:,2)
          q2=dot_product(v,v)
          ve(:,i)=v/sqrt(q2)
        enddo
        allocate(a(3,3))
        do i=1,3
          a(i,1)=ve(1,i)**2
          a(i,2)=2.d0*ve(1,i)*ve(2,i)
          a(i,3)=ve(2,i)**2
        enddo
        deallocate(ve)
        allocate(ipiv(3))
        call dgesv(3,ndim3_nu,a,3,ipiv,ph0,4,info)
        deallocate(a,ipiv)
        do ind_nu=1,ndim3_nu
          i=ind_tens33(1,1)
          pha(i,ind_nu)=ph0(1,ind_nu)
          i=ind_tens33(1,2)
          pha(i,ind_nu)=ph0(2,ind_nu)
          i=ind_tens33(2,1)
          pha(i,ind_nu)=ph0(2,ind_nu)
          i=ind_tens33(2,2)
          pha(i,ind_nu)=ph0(3,ind_nu)
          i=ind_tens33(3,3)
          pha(i,ind_nu)=ph0(4,ind_nu)
        enddo
      else if(istruc==121) then
        allocate(ve(3,3))
        do i=2,4
          if(i==2) v=qb0(:,2)
          if(i==3) v=qb0(:,2)+qb0(:,3)
          if(i==4) v=qb0(:,3)
          q2=dot_product(v,v)
          ve(:,i-1)=v/sqrt(q2)
        enddo
        allocate(a(3,3))
        do i=1,3
          a(i,1)=ve(1,i)**2
          a(i,2)=2.d0*ve(1,i)*ve(2,i)
          a(i,3)=ve(2,i)**2
        enddo
        deallocate(ve)
        allocate(ipiv(3))
        call dgesv(3,ndim3_nu,a,3,ipiv,ph0(2,1),4,info)
        deallocate(a,ipiv)
        do ind_nu=1,ndim3_nu
          i=ind_tens33(1,1)
          pha(i,ind_nu)=ph0(1,ind_nu)
          i=ind_tens33(2,2)
          pha(i,ind_nu)=ph0(2,ind_nu)
          i=ind_tens33(2,3)
          pha(i,ind_nu)=ph0(3,ind_nu)
          i=ind_tens33(3,2)
          pha(i,ind_nu)=ph0(3,ind_nu)
          i=ind_tens33(3,3)
          pha(i,ind_nu)=ph0(4,ind_nu)
        enddo
      else
        allocate(ve(3,6))
        do i=1,6
          if(i==1) v=qb0(:,1)
          if(i==2) v=qb0(:,2)
          if(i==3) v=qb0(:,3)
          if(i==4) v=qb0(:,1)+qb0(:,2)
          if(i==5) v=qb0(:,1)+qb0(:,3)
          if(i==6) v=qb0(:,2)+qb0(:,3)
          q2=dot_product(v,v)
          ve(:,i)=v/sqrt(q2)
        enddo
        allocate(a(6,6))
        do i=1,6
          a(i,1)=ve(1,i)**2
          a(i,2)=2.d0*ve(1,i)*ve(2,i)
          a(i,3)=ve(2,i)**2
          a(i,4)=2.d0*ve(1,i)*ve(3,i)
          a(i,5)=2.d0*ve(2,i)*ve(3,i)
          a(i,6)=ve(3,i)**2
        enddo
        deallocate(ve)
        allocate(ipiv(6))
        call dgesv(6,ndim3_nu,a,6,ipiv,ph0,6,info)
        deallocate(a,ipiv)
        do ind_nu=1,ndim3_nu
          i=ind_tens33(1,1)
          pha(i,ind_nu)=ph0(1,ind_nu)
          i=ind_tens33(1,2)
          pha(i,ind_nu)=ph0(2,ind_nu)
          i=ind_tens33(2,2)
          pha(i,ind_nu)=ph0(3,ind_nu)
          i=ind_tens33(1,3)
          pha(i,ind_nu)=ph0(4,ind_nu)
          i=ind_tens33(2,3)
          pha(i,ind_nu)=ph0(5,ind_nu)
          i=ind_tens33(3,3)
          pha(i,ind_nu)=ph0(6,ind_nu)
        enddo
      endif
      deallocate(ph0)
      end
