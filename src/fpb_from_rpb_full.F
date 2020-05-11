      subroutine fpb_from_rpb_full(pr,pf,k,n0r,n0f,key)
c     key  =  0 pf is assumed to be zero
c     key  =  1 the result will be added to pf
      use atom_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k,n0r,n0f,key
      complex*16, intent(in) :: pr(n0r,n0r)
      complex*16, intent(inout) :: pf(n0f,n0f)
      integer :: lget,ie,ie2,km,km2,iatom,isort,i_f,i_l,nlmi,lmi,li,nf,
     &           nl
      complex*16, allocatable :: t(:,:)
      if(key==0) pf=(0.d0,0.d0)
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
              t(1:nl,km2)=t(1:nl,km2)+rpb_fpb(ie,ie2,li,isort)
     &                               *pr(1:nl,km)
            enddo
          enddo
        enddo
      enddo
      if(nplwgw_red(k)/=0) then
        call zgemm('n','c',nl,nplwgw(k),nplwgw_red(k),(1.d0,0.d0),
     &             pr(1,n_pbmt_red+1),n0r,rpb_fpb_int(1,1,k),nplw_gw,
     &             (0.d0,0.d0),t(1,n_pbmt+1),nl)
      endif
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
              pf(km2,1:nf)=pf(km2,1:nf)+rpb_fpb(ie,ie2,li,isort)
     &                                 *t(km,1:nf)
            enddo
          enddo
        enddo
      enddo
      if(nplwgw_red(k)/=0) then
        call zgemm('n','n',nplwgw(k),nf,nplwgw_red(k),(1.d0,0.d0),
     &             rpb_fpb_int(1,1,k),nplw_gw,t(n_pbmt_red+1,1),nl,
     &             (1.d0,0.d0),pf(n_pbmt+1,1),n0f)
      endif
      deallocate(t)
      end
