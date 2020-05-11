      subroutine trans_pb_pw
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer :: k,n_pbt,ind_k,ind_gg(1),nn1,iatom,isort,i0,icase,i1,
     &           ind0,ind1,l,m,lm,ie,ii0,ii1,k5,i,ia,ib,ic,ii,ind,info,
     &           j,jj
      real*8 :: v(3),phase,pi2
      complex*16 :: c1,c2,cc1,cc2
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: s(:,:)
	  if(allocated(pw_pb)) deallocate(pw_pb)
      if(ndim3_k(me_k+1)/=0) allocate(pw_pb(n_pbtot,
     &                                      ndim3_k(me_k+1)))
      if(allocated(pw_pb_tild)) deallocate(pw_pb_tild)
      if(ndim3_k(me_k+1)/=0) allocate(pw_pb_tild(n_pbtot,
     &                                           ndim3_k(me_k+1)))
      allocate(s(nplw_gw,nplw_gw),ipiv(nplw_gw))
      ind_gg=1
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
c - Memory distribution for matrix-parallelization (current K-point) --
        n_pbt=n_pbmt+nplwgw(k)
        call trans_pb_pw_0(k,pw_pb(1,ind_k),n_pbt,ind_gg(1))
        pw_pb_tild(1:n_pbt,ind_k)=pw_pb(1:n_pbt,ind_k)
        do j=1,nplwgw(k)
          jj=indpw_gw(j,k)
          do i=1,nplwgw(k)
            ii=indpw_gw(i,k)
            ia=igbs(1,jj)-igbs(1,ii)
            ib=igbs(2,jj)-igbs(2,ii)
            ic=igbs(3,jj)-igbs(3,ii)
            ind=indplw(ia,ib,ic)
            if(complex_ro) then
              s(i,j)=dcmplx(sovr(ind),sovi(ind))
            else
              s(i,j)=dcmplx(sovr(ind),0.d0)
            endif
          enddo
        enddo
        if(nplwgw(k)/=0) call zgesv(nplwgw(k),1,s,nplw_gw,ipiv,
     &                              pw_pb_tild(n_pbmt+1,ind_k),
     &                              n_pbtot,info)
      enddo   !! over ind_k
      deallocate(s,ipiv)
	  if(.not.complex_ro) then
        pi2=pi+pi
	    if(allocated(pw_pb_r)) deallocate(pw_pb_r)
        if(allocated(pw_pb_tild_r)) deallocate(pw_pb_tild_r)
        if(ndim3_k(me_k+1)/=0) then
          allocate(pw_pb_r(n_pbtot,ndim3_k(me_k+1)))
          allocate(pw_pb_tild_r(n_pbtot,ndim3_k(me_k+1)))
        endif
        nn1=n_pbmt+1
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
	      n_pbt=n_pbmt+nplwgw(k)
          pw_pb_r(nn1:n_pbt,ind_k)=real(pw_pb(nn1:n_pbt,ind_k))
          pw_pb_tild_r(nn1:n_pbt,ind_k)=
     &             real(pw_pb_tild(nn1:n_pbt,ind_k))
          do iatom=1,natom
            isort=is(iatom)
            i0=ip(iatom,inv_num)
            if(i0==iatom) then
              icase=1
            else if(i0<iatom) then
              icase=2
              i1=iatom
            else if(i0>iatom) then
              icase=2
              i1=i0
              i0=iatom
            endif
            ind0=iopb(i0)-1
            if(icase==2) ind1=iopb(i1)-1
            if(icase==1) then
              v=tau(:,iatom)
              phase=pi2*dot_product(pnt(:,k),v)
              c1=dcmplx(cos(phase),-sin(phase))
            else if(icase==2) then
              v=tau(:,i0)
              phase=pi2*dot_product(pnt(:,k),v)
              c1=dcmplx(cos(phase),-sin(phase))/sqrt(2.d0)
              v=tau(:,i1)
              phase=pi2*dot_product(pnt(:,k),v)
              c2=dcmplx(cos(phase),-sin(phase))/sqrt(2.d0)
            endif
            do l=0,lmpb(isort)
              if(icase==1) then
                cc1=(0.d0,-1.d0)**l
              else if(icase==2) then
                cc1=(0.d0,-1.d0)**l
                cc2=(0.d0,-1.d0)**(l+1)
              endif
              do m=-l,l
                lm=l*(l+1)+m+1
                do ie=1,ntle_pb(l,isort)
                  k5=ind_prod(ie,lm,isort)
                  ii0=ind0+k5
                  if(icase==2) ii1=ind1+k5
                  if(icase==1) then
                    pw_pb_r(ii0,ind_k)=cc1*c1*pw_pb(ii0,ind_k)
                    pw_pb_tild_r(ii0,ind_k)=pw_pb_r(ii0,ind_k)
                  else if(icase==2) then
                    pw_pb_r(ii0,ind_k)=cc1*(c1*pw_pb(ii0,ind_k)
     &                                     +c2*pw_pb(ii1,ind_k))
                    pw_pb_r(ii1,ind_k)=cc2*(c2*pw_pb(ii1,ind_k)
     &                                     -c1*pw_pb(ii0,ind_k))
                    pw_pb_tild_r(ii0,ind_k)=pw_pb_r(ii0,ind_k)
                    pw_pb_tild_r(ii1,ind_k)=pw_pb_r(ii1,ind_k)
                  endif
                enddo   !! over ie
              enddo  !! over m
            enddo   !! over l
          enddo   !! over iatom
        enddo  !! over ind_k
      endif
      end
