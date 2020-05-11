      subroutine set_lattices
	  use solid_mod
      implicit none
	  integer :: nkmx
	  real*8, allocatable :: work(:),dlat0(:,:),rlat0(:,:)
	  if(allocated(dtau)) deallocate(dtau)
	  if(allocated(ipdtau)) deallocate(ipdtau)
	  allocate(dtau(3,natom**2),ipdtau(natom,natom))
      call convx(tau,ndtau)
	  nkmx=90000
	  allocate(dlat0(3,nkmx),rlat0(3,nkmx),work(nkmx))
	  call lattc(dlat0,rlat0,nkmx,nkmx,work)
	  if(allocated(dlat)) deallocate(dlat)
	  if(allocated(rlat)) deallocate(rlat)
	  allocate(dlat(3,nkd),rlat(3,nkr))
	  dlat=dlat0(:,1:nkd)
	  rlat=rlat0(:,1:nkr)
	  deallocate(work,dlat0,rlat0)
      end
      
      subroutine convx(tu,ndt)
c     collects pairs with equivalent connecting vector.
c     ipdtau(i,j).gt.0 points to tau_i-tau_j, .lt.0 to tau_j-tau_i.
c     don't forget bloch phase if equivalence is by a
c     lattice vector not equal to zero !
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
	  integer :: iatom,jatom,k,idold,ndt,id
      real*8 :: tu(3,natom),d1(3),d2(3),xx(3),tol
      tol=1.d-10
      ndt=0
      do iatom=1,natom
        do jatom=1,natom
          d1=tu(:,iatom)-tu(:,jatom)
          d2=-d1
          call zone0_cart(d1,rbas,d1)
          call zone0_cart(d2,rbas,d2)
          idold=0
          do k=1,ndt
            xx=d2-dtau(:,k)
            call zone0_cart(xx,rbas,xx)
            if(xx(1)**2+xx(2)**2+xx(3)**2.lt.tol) idold=-k
            xx=d1-dtau(:,k)
            call zone0_cart(xx,rbas,xx)
            if(xx(1)**2+xx(2)**2+xx(3)**2.lt.tol) idold=k
          enddo
          if(idold.eq.0) then
            ndt=ndt+1
            dtau(:,ndt)=d1
            idold=ndt
          endif
          ipdtau(iatom,jatom)=idold
        enddo
	  enddo
c ------------------------------------
      if(maswrk) write(iun,100) ndt
  100 format(' number of non-equivalent connecting vectors is',i4)
      if(maswrk) then
        do id=1,ndt
          write(iun,340) id,dtau(1,id),dtau(2,id),dtau(3,id)
        enddo
  340   format(i5,3x,3f12.5)
      endif
      end
