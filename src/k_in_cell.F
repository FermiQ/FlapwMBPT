	subroutine k_in_cell(q,k)
	use models_mod
	use solid_mod
	implicit none
	real*8, intent(in) :: q(3)
	integer, intent(out) :: k
      integer :: ij,i1,i2,i3
      real*8 :: q1(3),v(3),qb1,qb2,qb3,vb1,vb2,vb3
      call zone1_cart(q,gbas,q1)
      qb1=dot_product(qb0(:,1),qb0(:,1))
      qb2=dot_product(qb0(:,2),qb0(:,2))
      qb3=dot_product(qb0(:,3),qb0(:,3))
      ij=0
      do i3=0,ndiv(3)-1
        do i2=0,ndiv(2)-1
          do i1=0,ndiv(1)-1
            ij=ij+1
            v=qb0(:,1)*i1+qb0(:,2)*i2+qb0(:,3)*i3
            v=q1-v
            vb1=dot_product(v,qb0(:,1))/qb1
            vb2=dot_product(v,qb0(:,2))/qb2
            vb3=dot_product(v,qb0(:,3))/qb3
            if(vb1>=0.d0.and.vb1<1.d0) then
              if(vb2>=0.d0.and.vb2<1.d0) then
                if(vb3>=0.d0.and.vb3<1.d0) then
                  k=ij
                endif
              endif
            endif
          enddo
        enddo
      enddo
	end