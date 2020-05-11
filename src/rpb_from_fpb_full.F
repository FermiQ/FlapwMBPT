      subroutine rpb_from_fpb_full(pr,pf,k,n0r,n0f)
      use atom_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k,n0r,n0f
      complex*16, intent(in) :: pf(n0f,n0f)
      complex*16, intent(out) :: pr(n0r,n0r)
      integer :: lget,ie,ie2,km,km2,iatom,isort,i_f,i_l,nlmi,lmi,li,nf,
     &           nl
      complex*16, allocatable :: t(:,:)
      nf=n_pbmt+nplwgw(k)
      nl=n_pbmt_red+nplwgw_red(k)
      allocate(t(nl,nf))
      t=(0.d0,0.d0)
      do iatom=1,natom
        isort=is(iatom)
        i_f=iopb(iatom)-1
        i_l=iopb_red(iatom)-1
        nlmi=(lmpb_red(isort)+1)**2
        do lmi=1,nlmi
          li=lget(lmi)
          do ie2=1,ntle_pb(li,isort)
            km2=i_f+ind_prod(ie2,lmi,isort)
            do ie=1,ntle_pb_red(li,isort)
              km=i_l+ind_prod_red(ie,lmi,isort)
              t(km,1:nf)=t(km,1:nf)+rpb_fpb(ie,ie2,li,isort)
     &                             *pf(km2,1:nf)
            enddo
          enddo
        enddo
      enddo
      if(nplwgw_red(k)/=0) then
        call zgemm('c','n',nplwgw_red(k),nf,nplwgw(k),(1.d0,0.d0),
     &             rpb_fpb_int(1,1,k),nplw_gw,pf(n_pbmt+1,1),n0f,
     &             (0.d0,0.d0),t(n_pbmt_red+1,1),nl)
      endif
      pr=(0.d0,0.d0)
      do iatom=1,natom
        isort=is(iatom)
        i_f=iopb(iatom)-1
        i_l=iopb_red(iatom)-1
        nlmi=(lmpb_red(isort)+1)**2
        do lmi=1,nlmi
          li=lget(lmi)
          do ie2=1,ntle_pb(li,isort)
            km2=i_f+ind_prod(ie2,lmi,isort)
            do ie=1,ntle_pb_red(li,isort)
              km=i_l+ind_prod_red(ie,lmi,isort)
              pr(1:nl,km)=pr(1:nl,km)+rpb_fpb(ie,ie2,li,isort)
     &                               *t(1:nl,km2)
            enddo
          enddo
        enddo
      enddo
      if(nplwgw_red(k)/=0) then
        call zgemm('n','n',nl,nplwgw_red(k),nplwgw(k),(1.d0,0.d0),
     &             t(1,n_pbmt+1),nl,rpb_fpb_int(1,1,k),nplw_gw,
     &             (0.d0,0.d0),pr(1,n_pbmt_red+1),n0r)
      endif
      deallocate(t)
      end
