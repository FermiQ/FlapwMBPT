      subroutine sym_me_pb_red(fxc_mt)
      use atom_mod
      use manager_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer :: isort,iatom,ja,ig,jatom,j,l,lp1,ie,m,lm,km,lm0,i,i5
      real*8, allocatable :: y0(:),t0(:),pa1(:,:)
      real*8, intent(out) :: fxc_mt(n_pbmtm_red,n_pbmtm_red,natom,
     &                              n_ixc_0)
      allocate(y0(2*maxpb_red+1),t0(2*maxpb_red+1),
     &         pa1(n_pbmtm_red,n_pbmtm_red))
c ------ start loop over atoms in class ----------
      do isort=1,nsort
        iatom=iat_1(isort)
        do ja=iatom+1,natom     !!  over atoms of the given class
          do ig=1,ngroup
            jatom=ip(ja,ig)
            if(jatom.eq.iatom) then
              do i5=1,n_ixc_0
                do j=1,n_pbmt0_red(isort)
                  do l=0,lmpb_red(isort)
                    lp1=l+1
                    do ie=1,ntle_pb_red(l,isort)
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=ind_prod_red(ie,lm,isort)
                        lm0=l+m+1
                        t0(lm0)=fxc_mt(km,j,iatom,i5)
                      enddo
                      call rotate1(t0,y0,l,u(1,ig),1)
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=ind_prod_red(ie,lm,isort)
                        lm0=l+m+1
                        pa1(km,j)=y0(lm0)
                      enddo
                    enddo
                  enddo   !! over l
                enddo
                do i=1,n_pbmt0_red(isort)
                  do l=0,lmpb_red(isort)
                    lp1=l+1
                    do ie=1,ntle_pb_red(l,isort)
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=ind_prod_red(ie,lm,isort)
                        lm0=l+m+1
                        t0(lm0)=pa1(i,km)
                      enddo
                      call rotate1(t0,y0,l,u(1,ig),1)
                      do m=-l,l
                        lm=l*(l+1)+m+1
                        km=ind_prod_red(ie,lm,isort)
                        lm0=l+m+1
                        fxc_mt(i,km,ja,i5)=y0(lm0)
                      enddo
                    enddo
                  enddo
                enddo  !! over i
              enddo   !! over i5
              goto 2
            endif
          enddo   !! over ig
2         continue
        enddo   !! over ja
      enddo  !! over isort
      deallocate(y0,t0,pa1)
      end
