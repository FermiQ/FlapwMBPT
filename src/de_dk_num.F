      subroutine de_dk_num
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: i,ka,ka1,kb,kb1,kc,kc1,ispin,k
      real*8 :: qb(3,3),dd,fa,fb,fc,v(3),pi2,pi2a,det3
      if(.not.allocated(de_bnd)) allocate(de_bnd(nbndf,npnt,nspin_0,3))
      pi2=pi+pi
      pi2a=pi2/par
      qb=qb0*pi2a
      dd=det3(qb(1,1),qb(1,2),qb(1,3))
      do ispin=1,nspin
        do k=1,npnt
          v=pnt(:,k)+qb0(:,1)
          call zone1_number(v,rb0,ndiv,ka1)
          ka1=index_k1(ka1)
          ka1=i_kref(ka1)
          v=pnt(:,k)-qb0(:,1)
          call zone1_number(v,rb0,ndiv,ka)
          ka=index_k1(ka)
          ka=i_kref(ka)
          v=pnt(:,k)+qb0(:,2)
          call zone1_number(v,rb0,ndiv,kb1)
          kb1=index_k1(kb1)
          kb1=i_kref(kb1)
          v=pnt(:,k)-qb0(:,2)
          call zone1_number(v,rb0,ndiv,kb)
          kb=index_k1(kb)
          kb=i_kref(kb)
          v=pnt(:,k)+qb0(:,3)
          call zone1_number(v,rb0,ndiv,kc1)
          kc1=index_k1(kc1)
          kc1=i_kref(kc1)
          v=pnt(:,k)-qb0(:,3)
          call zone1_number(v,rb0,ndiv,kc)
          kc=index_k1(kc)
          kc=i_kref(kc)
          do i=1,n_bnd(k,ispin)
            fa=(e_bnd(i,ka1,ispin)-e_bnd(i,ka,ispin))/2
            fb=(e_bnd(i,kb1,ispin)-e_bnd(i,kb,ispin))/2
            fc=(e_bnd(i,kc1,ispin)-e_bnd(i,kc,ispin))/2
            de_bnd(i,k,ispin,1)=(fa*qb(2,2)*qb(3,3)+fc*qb(2,1)*qb(3,2)
     &                          +fb*qb(3,1)*qb(2,3)-fc*qb(3,1)*qb(2,2)
     &                       -fa*qb(3,2)*qb(2,3)-fb*qb(2,1)*qb(3,3))/dd
            de_bnd(i,k,ispin,2)=(fb*qb(1,1)*qb(3,3)+fa*qb(1,3)*qb(3,2)
     &                          +fc*qb(1,2)*qb(3,1)-fb*qb(1,3)*qb(3,1)
     &                       -fa*qb(1,2)*qb(3,3)-fc*qb(3,2)*qb(1,1))/dd
            de_bnd(i,k,ispin,3)=(fc*qb(2,2)*qb(1,1)+fa*qb(1,2)*qb(2,3)
     &                          +fb*qb(1,3)*qb(2,1)-fa*qb(1,3)*qb(2,2)
     &                       -fb*qb(1,1)*qb(2,3)-fc*qb(2,1)*qb(1,2))/dd
          enddo
        enddo
      enddo
      end


      subroutine de_dk_num_0(eb,deb,n)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: eb(n,npnt,nspin)
      real*8, intent(out) :: deb(n,npnt,nspin,3)
      integer :: i,ka,ka1,kb,kb1,kc,kc1,ispin,k
      real*8 :: qb(3,3),dd,fa,fb,fc,v(3),pi2,pi2a,det3
      pi2=pi+pi
      pi2a=pi2/par
      qb=qb0*pi2a
      dd=det3(qb(1,1),qb(1,2),qb(1,3))
      do ispin=1,nspin
        do k=1,npnt
          v=pnt(:,k)+qb0(:,1)
          call zone1_number(v,rb0,ndiv,ka1)
          ka1=index_k1(ka1)
          ka1=i_kref(ka1)
          v=pnt(:,k)-qb0(:,1)
          call zone1_number(v,rb0,ndiv,ka)
          ka=index_k1(ka)
          ka=i_kref(ka)
          v=pnt(:,k)+qb0(:,2)
          call zone1_number(v,rb0,ndiv,kb1)
          kb1=index_k1(kb1)
          kb1=i_kref(kb1)
          v=pnt(:,k)-qb0(:,2)
          call zone1_number(v,rb0,ndiv,kb)
          kb=index_k1(kb)
          kb=i_kref(kb)
          v=pnt(:,k)+qb0(:,3)
          call zone1_number(v,rb0,ndiv,kc1)
          kc1=index_k1(kc1)
          kc1=i_kref(kc1)
          v=pnt(:,k)-qb0(:,3)
          call zone1_number(v,rb0,ndiv,kc)
          kc=index_k1(kc)
          kc=i_kref(kc)
          do i=1,n
            fa=(eb(i,ka1,ispin)-eb(i,ka,ispin))/2
            fb=(eb(i,kb1,ispin)-eb(i,kb,ispin))/2
            fc=(eb(i,kc1,ispin)-eb(i,kc,ispin))/2
            deb(i,k,ispin,1)=(fa*qb(2,2)*qb(3,3)+fc*qb(2,1)*qb(3,2)
     &                       +fb*qb(3,1)*qb(2,3)-fc*qb(3,1)*qb(2,2)
     &                       -fa*qb(3,2)*qb(2,3)-fb*qb(2,1)*qb(3,3))/dd
            deb(i,k,ispin,2)=(fb*qb(1,1)*qb(3,3)+fa*qb(1,3)*qb(3,2)
     &                       +fc*qb(1,2)*qb(3,1)-fb*qb(1,3)*qb(3,1)
     &                       -fa*qb(1,2)*qb(3,3)-fc*qb(3,2)*qb(1,1))/dd
            deb(i,k,ispin,3)=(fc*qb(2,2)*qb(1,1)+fa*qb(1,2)*qb(2,3)
     &                       +fb*qb(1,3)*qb(2,1)-fa*qb(1,3)*qb(2,2)
     &                       -fb*qb(1,1)*qb(2,3)-fc*qb(2,1)*qb(1,2))/dd
          enddo
        enddo
      enddo
      end
