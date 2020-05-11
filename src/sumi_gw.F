      subroutine sumi_gw(correlation,ind_k,k,ispin,rhoi)
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      logical, intent(in) :: correlation
      integer, intent(in) :: ind_k,k,ispin
      complex*16, intent(inout) :: rhoi(nplwro,nspin)
      integer :: ibas,jbas,ind,jnd,ia,ib,ic,j,igb,ii,jgb,ipp,ih,
     &	         jp,jh,i,ndim,nbas,ibas0,irl,n0,jbas0
      real*8 :: cr(3),gki(3),gkj(3),wei,pi2,pi2a,c05,c2,g2i,gi,eni,
     &	        c2i,const1,anormi,g2j,gj,enj,c2j,anormj,const,const2,
     &          scal,ddot,ann,sr,si,womr(3),womi(3),relr,reli,de,green0,
     &          ee
      complex*16 :: cc
      real*8, allocatable :: gf00(:)
	  complex*16, allocatable :: tmp(:,:),aga(:,:),gf0(:,:)
	  n0=nbndf
	  ndim=n_bnd(k,ispin)
      if(.not.correlation) then
        allocate(gf00(ndim))
        if(k_integral=='FD') then
          do i=1,ndim
            de=e_bnd(i,k,ispin)-chem_pot
            gf00(i)=-green0(de,betta_t)
          enddo
        else if(k_integral=='O1') then
          gf00=0.d0
          do ii=1,nn_k_int
            do i=1,ndim
              ee=e_bnd(i,k,ispin)
     &          +de_bnd(i,k,ispin,1)*k_int_vec(1,ii)
     &          +de_bnd(i,k,ispin,2)*k_int_vec(2,ii)
     &          +de_bnd(i,k,ispin,3)*k_int_vec(3,ii)
              de=ee-chem_pot
              gf00(i)=gf00(i)-green0(de,betta_t)
            enddo
          enddo
          gf00=gf00/nn_k_int
        endif
      else
        allocate(gf0(ndim,ndim))
        gf0=(0.d0,0.d0)
        if(me_t==0) then
          call ferm_unpack_tau(gf0,g_full(1,1,1,1,ind_k,ispin),ndim,
     &                         ndim,nbndf,2)
          gf0=-gf0
        endif
        if(nproc_t/=1) call dgop(gf0,2*ndim**2,'  +',comm_t)
      endif
	  wei=wgt(k)*2.d0/dfloat(nrel*nspin)
      if(irel<=1) then
        nbas=nbask(k)
	    allocate(tmp(nbas,ndim),aga(nbas,nbas))
	    if(.not.correlation) then
	      do j=1,ndim
	        do i=1,nbas
              tmp(i,j)=ev_bnd(i,j,ind_k,ispin)*gf00(j)
	        enddo
	      enddo
	    else
	      call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &	             ev_bnd(1,1,ind_k,ispin),nbasmpw,gf0,ndim,
     &               (0.d0,0.d0),tmp,nbas)
        endif
	    call zgemm('n','c',nbas,nbas,ndim,(1.d0,0.d0),tmp,nbas,
     &	           ev_bnd(1,1,ind_k,ispin),nbasmpw,(0.d0,0.d0),aga,
     &             nbas)
        deallocate(tmp)
        aga=aga/amega
c --------------------------------------------------------------------
        do ibas=1,nbas
          ind=indgb(ibas,k)
          do jbas=1,nbas
            jnd=indgb(jbas,k)
            ia=igbs(1,ind)-igbs(1,jnd)
            ib=igbs(2,ind)-igbs(2,jnd)
            ic=igbs(3,ind)-igbs(3,jnd)
            j=indplw(ia,ib,ic)
            if(j>nplwro) cycle
            rhoi(j,ispin)=rhoi(j,ispin)+wei*aga(ibas,jbas)
          enddo   !!! over jbas
        enddo   !!! over ibas
        deallocate(aga)
      else if(irel==2) then
        if(.not.rel_interst) then
          nbas=nbask(k)/2
	      allocate(tmp(nbas,ndim),aga(nbas,nbas))
	      do irl=1,2
	        ibas0=0
	        if(irl==2) ibas0=nbas
	        if(.not.correlation) then
	          do j=1,ndim
	            do i=1,nbas
	              tmp(i,j)=ev_bnd(ibas0+i,j,ind_k,ispin)*gf00(j)
	            enddo
	          enddo
	        else
	          call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &	                 ev_bnd(ibas0+1,1,ind_k,1),nbasmpw,gf0,ndim,
     &                   (0.d0,0.d0),tmp,nbas)
            endif
	        call zgemm('n','c',nbas,nbas,ndim,(1.d0,0.d0),tmp,nbas,
     &	               ev_bnd(ibas0+1,1,ind_k,1),nbasmpw,(0.d0,0.d0),
     &                 aga,nbas)
            aga=aga/amega
c --------------------------------------------------------------------
            do ibas=1,nbas
              ind=indgb(ibas,k)
              do jbas=1,nbas
                jnd=indgb(jbas,k)
                ia=igbs(1,ind)-igbs(1,jnd)
                ib=igbs(2,ind)-igbs(2,jnd)
                ic=igbs(3,ind)-igbs(3,jnd)
                j=indplw(ia,ib,ic)
                if(j>nplwro) cycle
                cc=wei*aga(ibas,jbas)
                rhoi(j,1)=rhoi(j,1)+cc
                if(magn.eq.2) then
                  spintr(3,j)=spintr(3,j)+(-1)**irl*cc
                endif
              enddo   !!! over jbas
            enddo   !!! over ibas
          enddo !! over irl
          if(magn.eq.2) then
	        ibas0=0
	        jbas0=nbas
	        if(.not.correlation) then
	          do j=1,ndim
	            do i=1,nbas
	              tmp(i,j)=ev_bnd(ibas0+i,j,ind_k,ispin)*gf00(j)
	            enddo
	          enddo
	        else
	          call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &	                 ev_bnd(ibas0+1,1,ind_k,1),nbasmpw,gf0,ndim,
     &                   (0.d0,0.d0),tmp,nbas)
            endif
	        call zgemm('n','c',nbas,nbas,ndim,(1.d0,0.d0),tmp,nbas,
     &	               ev_bnd(jbas0+1,1,ind_k,1),nbasmpw,(0.d0,0.d0),
     &                 aga,nbas)
            aga=aga/amega
c --------------------------------------------------------------------
            do ibas=1,nbas
              ind=indgb(ibas,k)
              do jbas=1,nbas
                jnd=indgb(jbas,k)
                ia=igbs(1,ind)-igbs(1,jnd)
                ib=igbs(2,ind)-igbs(2,jnd)
                ic=igbs(3,ind)-igbs(3,jnd)
                j=indplw(ia,ib,ic)
                if(j>nplwro) cycle
                cc=wei*aga(ibas,jbas)
                spintr(1,j)=spintr(1,j)+cc
                spintr(2,j)=spintr(2,j)-(0.d0,1.d0)*cc
                cc=wei*conjg(aga(jbas,ibas))
                spintr(1,j)=spintr(1,j)+cc
                spintr(2,j)=spintr(2,j)+(0.d0,1.d0)*cc
              enddo   !!! over jbas
            enddo   !!! over ibas
          endif
          deallocate(aga,tmp)
        else if(rel_interst) then
          nbas=nbask(k)
	      allocate(tmp(nbas,ndim),aga(nbas,nbas))
	      if(.not.correlation) then
	        do j=1,ndim
	          do i=1,nbas
	            tmp(i,j)=ev_bnd(i,j,ind_k,ispin)*gf00(j)
	          enddo
	        enddo
	      else
	        call zgemm('n','n',nbas,ndim,ndim,(1.d0,0.d0),
     &	               ev_bnd(1,1,ind_k,1),nbasmpw,gf0,ndim,
     &                 (0.d0,0.d0),tmp,nbas)
          endif
	      call zgemm('n','c',nbas,nbas,ndim,(1.d0,0.d0),tmp,nbas,
     &	             ev_bnd(1,1,ind_k,ispin),nbasmpw,(0.d0,0.d0),aga,
     &               nbas)
          deallocate(tmp)
          aga=aga/amega
          pi2=2.d0*pi
          pi2a=pi2/par
          c05=0.5d0*clight
          c2=clight*clight
          do igb=1,nbask(k)/2
            ind=indgb(igb,k)
            g2i=0.d0
            do ii=1,3
              gki(ii)=(pnt(ii,k)+gbs(ii,ind))*pi2a
              g2i=g2i+gki(ii)**2
            enddo
            gi=sqrt(g2i)
            eni=c05*(sqrt(c2+(2.d0*gi)**2)-clight)
            c2i=c2+eni
            const1=clight/c2i
            anormi=c2i/sqrt(c2i*c2i+(clight*gi)**2)
            do jgb=1,nbask(k)/2
              jnd=indgb(jgb,k)
              ia=igbs(1,ind)-igbs(1,jnd)
              ib=igbs(2,ind)-igbs(2,jnd)
              ic=igbs(3,ind)-igbs(3,jnd)
              j=indplw(ia,ib,ic)
              if(j>nplwro) cycle
              g2j=0.d0
              do ii=1,3
                gkj(ii)=(pnt(ii,k)+gbs(ii,jnd))*pi2a
                g2j=g2j+gkj(ii)**2
              enddo
              gj=sqrt(g2j)
              enj=c05*(sqrt(c2+(2.d0*gj)**2)-clight)
              c2j=c2+enj
              anormj=c2j/dsqrt(c2j*c2j+(clight*gj)**2)
              const=c2/c2i/c2j
              const2=clight/c2j
              scal=ddot(3,gki,1,gkj,1)
              call cross(cr,gki,gkj)
              ann=anormi*anormj
              do ipp=1,2
                ih=(ipp-1)*nbask(k)/2+igb
                do jp=1,2
                  jh=(jp-1)*nbask(k)/2+jgb
                  sr=dreal(aga(ih,jh))
                  si=dimag(aga(ih,jh))
                  call setrel(cr,scal,const,ipp,jp,ann,gki,gkj,womr,
     &                        womi,relr,reli,magn)
                  rhoi(j,1)=rhoi(j,1)+wei*(sr*relr-si*reli)
                  if(complex_ro) rhoi(j,1)=rhoi(j,1)+(0.d0,1.d0)*wei
     &                                              *(sr*reli+si*relr)
                  if(magn.eq.2) then
                    do i=1,3
                      spintr(i,j)=spintr(i,j)+wei
     &                 *dcmplx(sr*womr(i)-si*womi(i),
     &                         sr*womi(i)+si*womr(i))
                    enddo
                  endif
                enddo  !!! jp
              enddo  !!! ipp
            enddo   !!! over jgb
          enddo   !!! over igb
          deallocate(aga)
        endif
      endif   !!! over irel
      if(correlation) deallocate(gf0)
      if(.not.correlation) deallocate(gf00)
      end
