      subroutine sbzmesh(kpnt,kwgt,nk,ierror)
	  use manager_mod
	  use parallel_mod
      use sitegen_mod
	  use solid_mod
      implicit none
	  integer, intent(out) :: nk,ierror
	  real*8, intent(out) :: kpnt(3,nqdiv),kwgt(nqdiv)
      integer :: nvcmax,n1,n2,n3,iz,jz,kz,i,nvec,ivec,
     &	         nwgt,i1,i2,i3,j1,j2,j3,k1,k2,k3,nred,iwgt,ig
      real*8 :: v(3),v1(3),v2(3),rmax,rsum,x1,x2,x3,scal,scal1,
     &          d_err,tu(3)
      integer, allocatable :: igq(:,:,:),ipq(:,:,:),iref(:),kvec(:,:)
      real*8, allocatable :: gvec(:,:)
      tu=0.d0
      ierror=0
      nvcmax=30000
      allocate(lastg_site(nvcmax/2))
      allocate(iref_vgen(nvcmax))
      allocate(vect_site(3,0:nvcmax))
      allocate(rco_site(nvcmax/2))
      allocate(tauc_site(3,natom))
      allocate(vec_vgen(3,nvcmax))
      allocate(rw_vgen(nvcmax))
c  set divisions along x,y,z
      n1=ndiv(1)-1
      n2=ndiv(2)-1
      n3=ndiv(3)-1
      allocate(igq(0:n1,0:n2,0:n3))
      allocate(ipq(0:n1,0:n2,0:n3))
      allocate(iref(0:nvcmax))
      allocate(kvec(3,0:nvcmax))
      allocate(gvec(3,0:nvcmax))
c  find sphere radius circumscribing bz
      rmax=0.d0
      do iz=0,1
        do jz=0,1
          do kz=0,1
            rsum=0.d0
            do i=1,3
              rsum=rsum+(iz*gbas(i,1)+jz*gbas(i,2)+kz*gbas(i,3))**2
            enddo
            rmax=max(rmax,rsum)
          enddo
        enddo
      enddo
      rmax=1.5d0*sqrt(rmax)+0.1d0
c  generate all reciprocal lattice vectors inside sphere "rmax"
      call sitegen(rmax,gbas(1,1),gbas(1,2),gbas(1,3),tu,1,nvec,gvec,
     &             iref,nvcmax)
      nvec=nvec-1
      do ivec=0,nvec
        x1=gvec(1,ivec)*rb0(1,1)+gvec(2,ivec)*rb0(2,1)+
     &                          gvec(3,ivec)*rb0(3,1)
        x2=gvec(1,ivec)*rb0(1,2)+gvec(2,ivec)*rb0(2,2)+
     &                          gvec(3,ivec)*rb0(3,2)
        x3=gvec(1,ivec)*rb0(1,3)+gvec(2,ivec)*rb0(2,3)+
     &                          gvec(3,ivec)*rb0(3,3)
        kvec(1,ivec)=anint(x1)
        kvec(2,ivec)=anint(x2)
        kvec(3,ivec)=anint(x3)
        if(abs(kvec(1,ivec)-x1).ge.1.d-05.or.
     &     abs(kvec(2,ivec)-x2).ge.1.d-05.or.
     &     abs(kvec(3,ivec)-x3).ge.1.d-05) then
          if(maswrk) write(iun,*)'<bzmesh>:
     &  error in reciprocal lattice vectors'
          call ending
        endif
      enddo
c  initialize point-indicator array
      igq=0
      ipq=0
c  *** find unequivalent k-points in whole bz ***
      nk=0
      nwgt=0
c  do over whole brillouin zone
      do i3=0,n3
        do i2=0,n2
          do i1=0,n1
            if(ipq(i1,i2,i3).eq.0)then
              v(1)=i1*qb0(1,1)+i2*qb0(1,2)+i3*qb0(1,3)
              v(2)=i1*qb0(2,1)+i2*qb0(2,2)+i3*qb0(2,3)
              v(3)=i1*qb0(3,1)+i2*qb0(3,2)+i3*qb0(3,3)
c     *************************************
              scal=dot_product(v,v)
              v2=v
c     *************************************
              iwgt=0
c  rotate given point
              do ig=1,ngroup
                call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),2)
                d_err=v1(1)*rb0(1,1)+v1(2)*rb0(2,1)+v1(3)*rb0(3,1)
                j1=anint(d_err)
                if(abs(dfloat(j1)-d_err)>1.d-6) ierror=1
                d_err=v1(1)*rb0(1,2)+v1(2)*rb0(2,2)+v1(3)*rb0(3,2)
                j2=anint(d_err)
                if(abs(dfloat(j2)-d_err)>1.d-6) ierror=2
                d_err=v1(1)*rb0(1,3)+v1(2)*rb0(2,3)+v1(3)*rb0(3,3)
                j3=anint(d_err)
                if(abs(dfloat(j3)-d_err)>1.d-6) ierror=3
c  reduce the point to the brillouin zone
                nred=0
                do ivec=0,nvec
                  k1=j1+kvec(1,ivec)
                  k2=j2+kvec(2,ivec)
                  k3=j3+kvec(3,ivec)
                  if((k1>=0.and.k1<=n1).and.(k2>=0.and.k2<=n2).and.
     &                                 (k3>=0.and.k3<=n3)) then
c  the point is reduced to the brillouin zone by vector ivec
                    nred=nred+1
                    if(nred.gt.8) then
                      if(maswrk) write(iun,*)
     &                   '<bzmesh>: # of reduced points is too large'
                      call ending
                    endif
                    if(igq(k1,k2,k3).eq.0) igq(k1,k2,k3)=ig
                    if(ipq(k1,k2,k3).eq.0) iwgt=iwgt+1
                    ipq(k1,k2,k3)=nk+1
                  endif
c     *************************************
                  v1(1)=k1*qb0(1,1)+k2*qb0(1,2)+k3*qb0(1,3)
                  v1(2)=k1*qb0(2,1)+k2*qb0(2,2)+k3*qb0(2,3)
                  v1(3)=k1*qb0(3,1)+k2*qb0(3,2)+k3*qb0(3,3)
                  scal1=dot_product(v1,v1)
                  if(scal1.lt.scal) then
                    scal=scal1
                    v2=v1
                  endif
c     *************************************
                enddo
                if(nred.eq.0)then
                  if(maswrk) write(iun,*)
     &            '<bzmesh>: the point is not reduced,
     &            check out for bz-translations and divisions n1,n2,n3'
                  call ending
                endif
              enddo
              nk=nk+1
              kpnt(1,nk)=v2(1)
              kpnt(2,nk)=v2(2)
              kpnt(3,nk)=v2(3)
              kwgt(nk)=iwgt
              nwgt=nwgt+iwgt
            endif
          enddo
        enddo
      enddo
      if(nwgt.ne.nqdiv)then
        if(maswrk) write(iun,*)'<tetra>: error in points'
        call ending
      endif
      kwgt=kwgt/nqdiv
	  deallocate(igq,ipq,iref,kvec,gvec)
      deallocate(lastg_site,vect_site,rco_site,tauc_site,iref_vgen,
     &           vec_vgen,rw_vgen)
      end
