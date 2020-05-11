      subroutine theta_mtm_red(fxc_mt)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: isort,iatom,nlmf,j0,jfn,i0,ifn,lm2,isym,mt,ii,
     &         irad,j,lmj,i,lmi,ix,icg1,icg2,icg
      real*8 :: cgg,dqdall,sqpi4
      real*8, intent(out) :: fxc_mt(n_pbmtm_red,n_pbmtm_red,natom,
     &                              n_ixc_0)
      real*8, allocatable :: t(:,:,:,:),work(:),work0(:)
      fxc_mt=0.d0
      allocate(t(maxlfpb_red,maxlfpb_red,maxnsym,n_ixc_0))
      allocate(work(0:maxnrad))
      allocate(work0(0:maxnrad))
      sqpi4=sqrt(4.d0*pi)
      do isort=1,nsort
        iatom=iat_1(isort)
        nlmf=(lmpb_red(isort)+1)**2
        do j0=1,lfun_pb_red(isort)
          jfn=ind_pb_red(j0,isort)
          do i0=1,lfun_pb_red(isort)
            ifn=ind_pb_red(i0,isort)
            do irad=0,nrad(isort)
              work0(irad)=fun_pb_red(irad+ifn)*fun_pb_red(irad+jfn)
     &                                        *dr(irad,isort)
     &                                        *r(irad,isort)**2
            enddo
            do lm2=1,nlmf
              if(sym(lm2,isort)) then
                isym=lmsym(lm2,isort)
                mt=indmt(isym,isort,1)
                do i=1,n_ixc_0
                  do irad=0,nrad(isort)
                    work(irad)=work0(irad)*theta_mt(mt+irad,i)
                  enddo
                  t(i0,j0,isym,i)=dqdall(h(isort),work,nrad(isort))
                enddo
              endif
            enddo   !! over lm2
          enddo   !! over i0
        enddo   !! over j0
        do j=1,n_pbmt0_red(isort)
          j0=lf_pb_red(j,isort)
          lmj=lm_pbmt_red(j,isort)
          do i=1,n_pbmt0_red(isort)
            i0=lf_pb_red(i,isort)
            lmi=lm_pbmt_red(i,isort)
            ix=max0(lmi,lmj)
            ix=(ix*(ix-1))/2+min0(lmi,lmj)
            icg1=indxcg(ix)
            icg2=indxcg(ix+1)-1
            do icg=icg1,icg2
              lm2 = jcg(icg)
              cgg=cg(icg)
              if(lm2.le.nlmf) then
                if(sym(lm2,isort)) then
                  isym=lmsym(lm2,isort)
                  do ii=1,n_ixc_0
                    fxc_mt(i,j,iatom,ii)=fxc_mt(i,j,iatom,ii)
     &                                   +cgg*t(i0,j0,isym,ii)
                  enddo
                endif
              endif
            enddo  !! over icg
          enddo   !! over i
        enddo   !! over j
      enddo  !! over isort
      deallocate(t,work,work0)
c ----------------------------------------------------------------
      call sym_me_pb_red(fxc_mt)
      end
