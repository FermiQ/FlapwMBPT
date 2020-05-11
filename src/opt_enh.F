      subroutine opt_enh(b11,b12,b21,b22,nne,nrax,thet,e)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use units_mod
	  use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nne,nrax
      complex*16, intent(in) :: b11(0:nrax,3,3),
     &                          b21(nne,0:nrax,3),b12(nne,0:nrax,3),
     &                          b22(nne,nne,0:nrax),thet(nne,nne)
      complex*16, intent(out) :: e(3,3,0:nrax)
      integer :: inu,it,j,i
      real*8 :: dir(3),sq2
      complex*16 :: bb11,x11,y11,r11,e11(6)
      complex*16, allocatable :: bb12(:),bb21(:),x12(:),x21(:),x22(:,:),
     &                           tmp(:,:),y12(:),y21(:),y22(:,:)
      sq2=sqrt(2.d0)/2.d0
      allocate(bb12(nne))
      allocate(bb21(nne))
      allocate(x12(nne))
      allocate(x21(nne))
      allocate(x22(nne,nne))
      allocate(tmp(nne,nne))
      allocate(y12(nne))
      allocate(y21(nne))
      allocate(y22(nne,nne))
      do inu=0,nrax
        do it=1,6
          if(it==1) dir=(/1.d0,0.d0,0.d0/)
          if(it==2) dir=(/0.d0,1.d0,0.d0/)
          if(it==3) dir=(/0.d0,0.d0,1.d0/)
          if(it==4) dir=(/sq2,sq2,0.d0/)
          if(it==5) dir=(/sq2,0.d0,sq2/)
          if(it==6) dir=(/0.d0,sq2,sq2/)
          bb11=(0.d0,0.d0)
          do j=1,3
            do i=1,3
              bb11=bb11+dir(i)*dir(j)*b11(inu,i,j)
            enddo
          enddo
          bb12=(0.d0,0.d0)
          bb21=(0.d0,0.d0)
          do i=1,3
            bb12=bb12+dir(i)*b12(:,inu,i)
            bb21=bb21+dir(i)*b21(:,inu,i)
          enddo
c --------------    1-P[V+Theta] -------------------------------          
          x11=1.d0-8.d0*pi*bb11
          call zgemm('n','n',1,nne,nne,(-1.d0,0.d0),bb12,1,thet,nne,
     &               (0.d0,0.d0),x12,1)
          x21=-8.d0*pi*bb21
          call zgemm('n','n',nne,nne,nne,(-1.d0,0.d0),b22(1,1,inu),nne,
     &               thet,nne,(0.d0,0.d0),x22,nne)
          do i=1,nne
            x22(i,i)=x22(i,i)+(1.d0,0.d0)
          enddo
c --------------    (1-P[V+Theta])^(-1) -------------------------------
          do j=1,nne
            do i=1,nne
              y22(i,j)=x22(i,j)-x21(i)*x12(j)/x11
            enddo
          enddo
          call invers_z(nne,y22,nne)
          call zgemm('n','n',1,nne,nne,(-1.d0,0.d0),x12,1,y22,nne,
     &               (0.d0,0.d0),y12,1)
          y12=y12/x11
          call invers_z(nne,x22,nne)
          call zgemm('n','n',nne,1,nne,(1.d0,0.d0),x22,nne,x21,nne,
     &               (0.d0,0.d0),y21,nne)
          y11=x11
          do i=1,nne
            y11=y11-x12(i)*y21(i)
          enddo
          y11=(1.d0,0.d0)/y11
          x21=x21*y11
          call zgemm('n','n',nne,1,nne,(-1.d0,0.d0),x22,nne,x21,nne,
     &               (0.d0,0.d0),y21,nne)
c -------------------  R11 -------------------------------------
          r11=y11*bb11
          do i=1,nne
            r11=r11+y12(i)*bb21(i)
          enddo
c ------------------------ E00^(-1) -------------------------------
          e11(it)=1.d0+8.d0*pi*r11
        enddo  !! over it
        e(1,1,inu)=(1.d0,0.d0)/e11(1)
        e(2,2,inu)=(1.d0,0.d0)/e11(2)
        e(3,3,inu)=(1.d0,0.d0)/e11(3)
        e(1,2,inu)=(1.d0,0.d0)/e11(4)-0.5d0*(e(1,1,inu)+e(2,2,inu))
        e(1,3,inu)=(1.d0,0.d0)/e11(5)-0.5d0*(e(1,1,inu)+e(3,3,inu))
        e(2,3,inu)=(1.d0,0.d0)/e11(6)-0.5d0*(e(2,2,inu)+e(3,3,inu))
        e(2,1,inu)=e(1,2,inu)
        e(3,1,inu)=e(1,3,inu)
        e(3,2,inu)=e(2,3,inu)
      enddo
      deallocate(bb12,bb21,x12,x21,x22,tmp,y12,y21,y22)
      end
