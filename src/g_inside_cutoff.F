      subroutine g_inside_cutoff(num,pt,cut,order,ind,indx,n,n0,nd)
c     num - max index of GBS vector to look ---
c     pt(3) - k-point
c     cut   - cutoff
c     order=T if we order the vectors (ascending)
c     ind=T if we need index array INDX(n0)      
	  use manager_mod
	  use solid_mod
	  use units_mod
      implicit none
      logical, intent(in) :: order,ind
      integer, intent(in) :: num,n0,nd(3)
      real*8, intent(in) :: pt(3),cut
      integer, intent(out) :: n,indx(n0)
	  integer :: ig,i0,i1,i,nr
	  real*8 :: ek,gk,p(3)
	  integer, allocatable :: flag(:),i00(:),i11(:)
	  real*8, allocatable :: gb(:,:),g0(:,:)
	  allocate(gb(2,n0))
	  allocate(flag(num))
	  n=0
c  ---- First include all vectors <= CUTOFF --------
      do ig=1,num
        p=pt+gbs(:,ig)
        ek=dot_product(p,p)
        gk=sqrt(ek)
        if(gk<=cut) then
          n=n+1
          if(order) then
            do i0=1,n-1
              if(gb(2,i0)>gk) then
                do i1=n,i0+1,-1
                  do i=1,2
                    gb(i,i1)=gb(i,i1-1)
                  enddo
                enddo
                gb(1,i0)=dfloat(ig)
                gb(2,i0)=gk
                goto 1
              endif
            enddo
          endif
          gb(1,n)=dfloat(ig)
          gb(2,n)=gk
        endif
1       continue
      enddo
c  ---- Now we remove periodically repeating images ---
      nr=nd(1)*nd(2)*nd(3)
      allocate(i00(n))
      allocate(i11(nr))
      i11=0
      do i=1,n
        ig=idnint(gb(1,i))
        call zone1_number(gbs(:,ig),rbas,nd,i00(i))
        i11(i00(i))=i11(i00(i))+1
      enddo
      flag=0
      do i=1,nr
        if(i11(i)<=1) cycle
        do i0=1,n
          if(i00(i0)==i) flag(i0)=1
        enddo
      enddo
      allocate(g0(2,n))
      i0=0
      do i=1,n
        if(flag(i)==1) cycle
        i0=i0+1
        g0(:,i0)=gb(:,i)
      enddo
      n=i0
c ------------------------------------------------------
      if(ind) indx(1:n)=idnint(g0(1,1:n))
	  deallocate(gb,flag,i00,i11,g0)
      end
