      subroutine geometry
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      character*1 :: ap
      integer :: ig,jatom,iatom,ia,ib,ic,nat,i,isort,n,ipiv(3),info
      real*8 :: x(3),x1(3),t(3),sd(3),a(3,3)
      integer, allocatable :: is0(:)
      real*8, allocatable :: t0(:,:)
      ap="'"
c   Account for orthorombicity parameters
      do i=1,ngroup
         shift(2,i)=shift(2,i)*b_a
         shift(3,i)=shift(3,i)*c_a
      enddo
      do i=1,3
         rbas(2,i)=rbas(2,i)*b_a
         rbas(3,i)=rbas(3,i)*c_a
      enddo
      do i=1,nsort
         tau(2,i)=tau(2,i)*b_a
         tau(3,i)=tau(3,i)*c_a
      enddo
      allocate(nhsort(nsort))
      nhsort=0
      nat=nsort*ngroup
      allocate(is0(nat))
      allocate(t0(3,nat))
      natom=0
      do isort=1,nsort
        do ig=1,ngroup
          call rotate(tau(1,isort),tau(2,isort),tau(3,isort),x(1),
     &                x(2),x(3),u(2,ig),2)
          x=x+shift(:,ig)
          call zone0_cart(x,rbas,x)
          n=0
          do iatom=1,natom
            x1=x-t0(:,iatom)
            do ia=-3,3
              do ib=-3,3
                do ic=-3,3
                  t(1)=ia*rbas(1,1)+ib*rbas(1,2)+ic*rbas(1,3)
                  t(2)=ia*rbas(2,1)+ib*rbas(2,2)+ic*rbas(2,3)
                  t(3)=ia*rbas(3,1)+ib*rbas(3,2)+ic*rbas(3,3)
                  sd=t-x1
                  if(abs(sd(1))+abs(sd(2))+abs(sd(3)).lt.1.d-6) n=1
                enddo
              enddo
            enddo
          enddo
          if(n==0) then
            natom=natom+1
            is0(iatom)=isort
            nhsort(isort)=nhsort(isort)+1
            t0(:,iatom)=x
          endif
        enddo
      enddo
      deallocate(tau)
      allocate(tau(3,natom))
      tau=t0(:,1:natom)
      allocate(is(natom))
      is=is0(1:natom)
      deallocate(is0,t0)
c ----- Check of symmetry operations relative to atomic positions --
      allocate(ip(natom,148))
      allocate(tshift(3,natom,ngroup))
      do ig=1,ngroup
        jat:do jatom=1,natom
          call rotate(tau(1,jatom),tau(2,jatom),tau(3,jatom),x(1),
     &                x(2),x(3),u(2,ig),2)   !!! = A*t^(-1)_a
          x=x+shift(:,ig)
          do iatom=1,natom
            x1=tau(:,iatom)-x
            do ia=-3,3
              do ib=-3,3
                do ic=-3,3
                  t(1)=ia*rbas(1,1)+ib*rbas(1,2)+ic*rbas(1,3)
                  t(2)=ia*rbas(2,1)+ib*rbas(2,2)+ic*rbas(2,3)
                  t(3)=ia*rbas(3,1)+ib*rbas(3,2)+ic*rbas(3,3)
                  sd=t-x1
                  if(abs(sd(1))+abs(sd(2))+abs(sd(3)).lt.1.d-6) then
                    ip(iatom,ig)=jatom
                    tshift(:,iatom,ig)=t
                    cycle jat
                  endif
                enddo
              enddo
            enddo
          enddo
          if(maswrk) then
            write(iun,*)'Not find jatom for ig =',ig,' iatom =',jatom
            write(iun,*)tau(1,jatom),' => ',x(1)
            write(iun,*)tau(2,jatom),' => ',x(2)
            write(iun,*)tau(3,jatom),' => ',x(3)
          endif
          call ending
        enddo jat
      enddo
      if(maswrk) then
        write(iun,*)' IS = '
        write(iun,'(18i4)')(is(iatom),iatom=1,natom)
        write(iun,*)
     &  ' Atomic Coordinates (in units of PAR, before ORT scaling)'
        do iatom=1,natom
          isort=is(iatom)
          write(iun,'(i5,1x,a4,3f20.15)')iatom,txtel(isort),
     &                   tau(1,iatom),tau(2,iatom)/b_a,tau(3,iatom)/c_a
        enddo
        write(iun,*)
     &  ' Atomic Coordinates (in units of PAR, after ORT scaling)'
        do iatom=1,natom
          isort=is(iatom)
          write(iun,'(i5,1x,a4,3f20.15)')iatom,txtel(isort),
     &                                   (tau(i,iatom),i=1,3)
        enddo
        write(iun,*)' Fractional Atomic Coordinates '
        do iatom=1,natom
          isort=is(iatom)
          x=tau(:,iatom)
          a=rbas
          call dgesv(3,1,a,3,ipiv,x,3,info)
          write(iun,'(i5,1x,a4,3f20.15)')iatom,txtel(isort),(x(i),i=1,3)
c          do ic=1,4
c            do ib=1,2
c              do ia=1,2
c                x1(1)=x(1)+ia-1
c                x1(2)=x(2)+ib-1
c                x1(3)=x(3)+ic-1
c                print '(3(a1,f5.2),a1)',ap,x1(1),',',x1(2),',',x1(3),ap
c              enddo
c            enddo
c          enddo
        enddo
      endif
c      stop
c ----- Check of symmetry operations relative to external m. field -
      if(irel==2.and.magn==2) call check_b
c ------------------------------------------------------------------
      nat_sort_max=maxval(nhsort)
      allocate(iat_sort(nat_sort_max,nsort))
      allocate(iat_sort_back(natom))
      do isort=1,nsort
        i=0
        do iatom=1,natom
          if(is(iatom)==isort) then
            i=i+1
            iat_sort(i,isort)=iatom
            iat_sort_back(iatom)=i
          endif
        enddo
      enddo
      allocate(iat_1(nsort))
      do isort=1,nsort
c --------- First atoms of the classes -----------------------------
        do iatom=1,natom
          if(is(iatom).eq.isort) then
            iat_1(isort)=iatom
            exit
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      par=par*v_v0**(1.d0/3.d0)
      call get_bz
      do isort=1,nsort
        if(smt(isort)>=0.d0) then
          smt(isort)=smt(isort)*v_v0**(1.d0/3.d0)
        else if(smt(isort)<0.d0) then
          smt(isort)=-smt(isort)
        endif
      enddo
      call timel('**** GEOMETRY finished *************')
      end
