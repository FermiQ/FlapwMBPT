      subroutine q_long_setup
      use atom_mod
      use manager_mod
      use solid_mod
      implicit none
      integer :: k,i,j,m,iq,iq0,ig,gbs_number
      real*8 :: v(3),v1(3),qstep(3),av,av1,d,ab,v2(3)
      integer, allocatable :: ns(:),ins(:,:,:)
      v=q1_long-q0_long
      av=sqrt(dot_product(v,v))
      d=1.d5
      do k=-5,5
        do j=-5,5
          do i=-5,5
            v1=i*qb0(:,1)+j*qb0(:,2)+k*qb0(:,3)
            av1=sqrt(dot_product(v1,v1))
            if(av1<1.d-6) cycle
            ab=dot_product(v,v1)
            if(abs(ab-av*av1)<1.d-6) then
              if(av1<d) then
                d=av1
                qstep=v1
              endif
            endif
          enddo
        enddo
      enddo
      step_q_long=d
      ab=av*q_long_factor
      m=idnint(ab/d)+1
      allocate(ns(npnt))
      allocate(ins(2,m,npnt))
      ns=0
      n_q_long=0
      do i=0,m
        if(d*i>ab) exit
        v=q0_long+qstep*i
        call zone1_number(v,rb0,ndiv,iq)
        iq=index_k1(iq)
        iq0=i_kref(iq)
        ns(iq0)=ns(iq0)+1
        n_q_long=n_q_long+1
        ins(2,ns(iq0),iq0)=n_q_long
        ig=k_group(iq)
        v1=v-pnt(:,iq)+gbs(:,g_sym_0(iq,ig))
	  call rotate(v1(1),v1(2),v1(3),v2(1),v2(2),v2(3),u(2,ig),1)
        ins(1,ns(iq0),iq0)=gbs_number(v2)
      enddo
      n_q0_long=0
      do k=1,npnt
        if(ns(k)/=0) n_q0_long=n_q0_long+1
      enddo
      maxsize_q_long=maxval(ns)
      allocate(q0_long_in_npnt(n_q0_long))
      allocate(size_q_long(n_q0_long))
      allocate(index_q_long(2,maxsize_q_long,n_q0_long))
      j=0
      do k=1,npnt
        if(ns(k)/=0) then
          j=j+1
          q0_long_in_npnt(j)=k
          size_q_long(j)=ns(k)
          do i=1,ns(k)
            index_q_long(:,i,j)=ins(:,i,k)
          enddo
        endif
      enddo
      deallocate(ns,ins)
      end