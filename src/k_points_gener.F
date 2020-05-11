      subroutine k_points_gener(nd,ndd,np,kp,wg)
      use manager_mod
      use parallel_mod
      use sitegen_mod
      use solid_mod
      implicit none
      integer, intent(in) :: nd(3),ndd
      integer, intent(out) :: np
      real*8, intent(out) :: kp(3,ndd),wg(ndd)
      integer :: n1,n2,n3,ind,ind1,i1,i2,i3,ig,m
      real*8 :: v(3),v1(3),rb(3,3),qb(3,3)
      integer, allocatable :: ind_s(:)
      do m = 1, 3
        rb(m,1) = rbas(m,1)*nd(1)
        rb(m,2) = rbas(m,2)*nd(2)
        rb(m,3) = rbas(m,3)*nd(3)
        qb(m,1) = gbas(m,1)/nd(1)
        qb(m,2) = gbas(m,2)/nd(2)
        qb(m,3) = gbas(m,3)/nd(3)
      enddo
      n1=nd(1)-1
      n2=nd(2)-1
      n3=nd(3)-1
      allocate(ind_s(ndd))
      ind_s=0
      np=0
      do i3=0,n3
        do i2=0,n2
          do i1=0,n1
            v=i1*qb(:,1)+i2*qb(:,2)+i3*qb(:,3)
            call zone1_number(v,rb,nd,ind)
            if(ind_s(ind)==0) then
              np=np+1
              kp(:,np)=v
              do ig=1,ngroup
                call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),2)
                call zone1_number(v1,rb,nd,ind1)
                ind_s(ind1)=np
              enddo
            endif
          enddo
        enddo
      enddo
      wg=0.d0
      do i1=1,ndd
        i2=ind_s(i1)
        wg(i2)=wg(i2)+1.d0
      enddo
      wg=wg/ndd
      end
