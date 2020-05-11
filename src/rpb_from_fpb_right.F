      subroutine rpb_from_fpb_right(pr,pf,nef,npf,npr,k)
      use atom_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nef,npf,npr,k
      complex*16, intent(in) :: pf(nef,npf)
      complex*16, intent(out) :: pr(nef,npr)
      integer :: lget,ie,ie2,km,km2,iatom,isort,i_f,i_l,nlmi,lmi,li
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
              pr(1:nef,km)=pr(1:nef,km)+rpb_fpb(ie,ie2,li,isort)
     &                                 *pf(1:nef,km2)
            enddo
          enddo
        enddo
      enddo
      if(nplwgw_red(k)/=0) then
        call zgemm('n','n',nef,nplwgw_red(k),nplwgw(k),(1.d0,0.d0),
     &             pf(1,n_pbmt+1),nef,rpb_fpb_int(1,1,k),nplw_gw,
     &             (0.d0,0.d0),pr(1,n_pbmt_red+1),nef)
      endif
      end
