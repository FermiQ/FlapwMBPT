      subroutine p_wing_num(pw,pw1)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: pw(ncmpl*n_pbtot,0:3,ndim3_nu)
      real*8, intent(out) :: pw1(ncmpl*3,n_opt_pb,ndim3_nu)
      integer :: i,ind_nu,iq,ipb,ipb0,n_pbt,npw
      real*8 :: pi2a,v(3),q1
      complex*16 :: cc(3),cc1(3)
      real*8, allocatable :: a(:,:),p1r(:,:),p0r(:,:)
      complex*16, allocatable :: p1(:,:),p0(:,:)
      pi2a=(pi+pi)/par
      npw=nplwgw(1)
      n_pbt=n_pbmt+npw
      if(ncmpl==1) then
        allocate(p1r(3,n_pbt))
      else
        allocate(p1(3,n_pbt))
      endif
      allocate(a(3,3))
      do i=1,3
        q1=sqrt(dot_product(qb0(:,i),qb0(:,i)))
        a(:,i)=qb0(:,i)/q1
      enddo
      do ind_nu=1,ndim3_nu
        if(ncmpl==1) then
          allocate(p0r(n_pbtot,0:3))
          call dcopy(n_pbtot*4,pw(1,0,ind_nu),1,p0r,1)
          p1r=0.d0
        else if(ncmpl==2) then
          allocate(p0(n_pbtot,0:3))
          call zcopy(n_pbtot*4,pw(1,0,ind_nu),1,p0,1)
          p1=(0.d0,0.d0)
        endif
        do i=1,3
          q1=pi2a*sqrt(dot_product(qb0(:,i),qb0(:,i)))
          call zone1_number(qb0(1,i),rb0,ndiv,iq)
          iq=index_k1(iq)
          iq=i_kref(iq)
          do ipb0=1,n_pbt
            if(ipb0<=n_pbmt) then
              ipb=ipb0
            else
              ipb=indpw_gw(ipb0-n_pbmt,1)
              ipb=iplf_gk(ipb,iq)
              if(ipb==0) goto 1
              ipb=n_pbmt+ipb
            endif
            if(ncmpl==1) p1r(i,ipb0)=(p0r(ipb,i)-p0r(ipb0,0))/q1
            if(ncmpl==2) p1(i,ipb0)=(p0(ipb,i)-p0(ipb0,0))/q1
1           continue
          enddo
        enddo
        if(ncmpl==1) deallocate(p0r)
        if(ncmpl==2) deallocate(p0)
        if(ncmpl==1) then
          call dgemm('n','t',3,n_opt_pb,n_pbt,1.d0,p1r,3,bm_coef_r,
     &               n_opt_pb,0.d0,pw1(1,1,ind_nu),3)
          do ipb=1,n_opt_pb
            v=pw1(:,ipb,ind_nu)
            do i=1,3
              pw1(:,ipb,ind_nu)=a(i,1)*v(1)+a(i,2)*v(2)+a(i,3)*v(3)
            enddo
          enddo
        else if(ncmpl==2) then
          call zgemm('n','c',3,n_opt_pb,n_pbt,(1.d0,0.d0),p1,3,bm_coef,
     &               n_opt_pb,(0.d0,0.d0),pw1(1,1,ind_nu),3)
          do ipb=1,n_opt_pb
            call zcopy(3,pw1(1,ipb,ind_nu),1,cc,1)
            do i=1,3
              cc1(i)=a(i,1)*cc(1)+a(i,2)*cc(2)+a(i,3)*cc(3)
            enddo
            call zcopy(3,cc1,1,pw1(1,ipb,ind_nu),1)
          enddo
        endif
      enddo   !! over ind_nu
      if(ncmpl==1) deallocate(p1r)
      if(ncmpl==2) deallocate(p1)
      deallocate(a)
      end
