      subroutine dos(dossp,pdos,weigd,ispin,ind_k,n)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      IMPLICIT none
      integer, intent(in) :: n,ispin,ind_k
      real*8, intent(in) :: weigd(n)
      real*8, intent(inout) :: pdos(nrel*maxb+1,nsort,nspin),
     &                         dossp(nspin)
      INTEGER :: k,ibnd,iatom,isort,ndimr,l2m2,li2,jj2,l2,
     &	         mj2,ie2,km2,ie1,km1,ind0,in,jn,in1,jn1,isp
      real*8, allocatable :: pdos0(:,:)
      allocate(pdos0(nrel*maxb+1,natom))
      isp=min(ispin,nspin_1)
      pdos0=0.d0
      k=n3_mpi_k(me_k+1)+ind_k
      do ibnd=1,n
        do iatom=1,natom
          ind0=io_lem(iatom)-1
          isort=is(iatom)
          ndimr=nrel*lmb(isort)+1
          l2m2=0
          do li2=1,ndimr
            if(irel.ne.2) then
              jj2=2*(li2-1)
              l2=li2-1
            else if(irel.eq.2) then
              jj2=li2+mod(li2,2)-1
              l2=li2/2
            endif
            do mj2=-jj2,jj2,2
              l2m2=l2m2+1
              do ie2=1,ntle(l2,isort)
                in1=1
                if(augm(ie2,l2,isort)=='LAW') in1=2
                if(augm(ie2,l2,isort)=='AWL') in1=2
                do jn1=1,in1
                  km2=ind0+indbasa(jn1,ie2,l2m2,isort)
                  do ie1=1,ntle(l2,isort)
                    in=1
                    if(augm(ie1,l2,isort)=='LAW') in=2
                    if(augm(ie1,l2,isort)=='AWL') in=2
                    do jn=1,in
                      km1=ind0+indbasa(jn,ie1,l2m2,isort)
                      pdos0(li2,iatom)=pdos0(li2,iatom)
     &		   		        +conjg(z_bnd(km1,ibnd,ind_k,ispin))
     &                            *z_bnd(km2,ibnd,ind_k,ispin)
     &                        *weigd(ibnd)
     &                        *ffsmt(jn,jn1,ie1,ie2,li2,isort,isp)
                    enddo
                  enddo   !!! over ie1
                enddo
              enddo   !!! over ie2
            enddo   !!! over mj2
          enddo   !!! over li2
        enddo        !!! over iatom
        dossp(ispin)=dossp(ispin)+weigd(ibnd)
      enddo               !!! over bands
      do isort=1,nsort
        ndimr=nrel*lmb(isort)+1
        do iatom=1,natom
          if(is(iatom)==isort) then
            do li2=1,ndimr
              pdos(li2,isort,ispin)=pdos(li2,isort,ispin)
     &	 	                       +pdos0(li2,iatom)
            enddo
          endif
        enddo
      enddo
      deallocate(pdos0)
      END
      
      
      
      
      subroutine pdos_x(pdos,pdos_at,pdos_diag,weigd,n0)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      IMPLICIT none
      integer, intent(in) :: n0
      real*8, intent(in) :: weigd(n0,npnt,nspin)
      real*8, intent(out) :: pdos(nrel*maxb+1,nsort,nspin),
     &                       pdos_at(0:nsort,nspin),
     &                       pdos_diag(maxntle,nrel*maxb+1,nsort,nspin) 
      INTEGER :: ispin,k,ind_k,ibnd,iatom,isort,ndimr,l2m2,li2,jj2,l2,
     &	       mj2,ie2,km2,ie1,km1,ndim,ind0,in,jn,in1,jn1,isp,l2m
      real*8 :: aa,cc,bb
      real*8, allocatable :: pdos0(:,:),pdos0_at(:),pdos0_diag(:,:,:)
      allocate(pdos0(nrel*maxb+1,natom))
      allocate(pdos0_at(0:natom))
      allocate(pdos0_diag(maxntle,nrel*maxb+1,natom))
      DO ISPIN=1,NSPIN
        isp=min(ispin,nspin_1)
        pdos0=0.d0
        pdos0_at=0.d0
        pdos0_diag=0.d0
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          ndim=n_bnd(k,ispin)
          do ibnd=1,ndim
            aa=0.d0
            do iatom=1,natom
              ind0=io_lem(iatom)-1
              isort=is(iatom)
              ndimr=nrel*lmb(isort)+1
              l2m=0
              do li2=1,ndimr
                if(irel.ne.2) then
                  jj2=2*(li2-1)
                  l2=li2-1
                else if(irel.eq.2) then
                  jj2=li2+mod(li2,2)-1
                  l2=li2/2
                endif
                do ie2=1,ntle(l2,isort)
                  in1=1
                  if(augm(ie2,l2,isort)=='LAW') in1=2
                  if(augm(ie2,l2,isort)=='AWL') in1=2
                  do jn1=1,in1
                    do ie1=1,ntle(l2,isort)
                      in=1
                      if(augm(ie1,l2,isort)=='LAW') in=2
                      if(augm(ie1,l2,isort)=='AWL') in=2
                      do jn=1,in
                        cc=0.d0
                        l2m2=l2m
                        do mj2=-jj2,jj2,2
                          l2m2=l2m2+1
                          km1=ind0+indbasa(jn,ie1,l2m2,isort)
                          km2=ind0+indbasa(jn1,ie2,l2m2,isort)
                          cc=cc+conjg(z_bnd(km1,ibnd,ind_k,ispin))
     &                               *z_bnd(km2,ibnd,ind_k,ispin)
                        enddo
                        bb=cc*ffsmt(jn,jn1,ie1,ie2,li2,isort,isp)
                        aa=aa+bb
                        pdos0(li2,iatom)=pdos0(li2,iatom)+bb
     &                        *weigd(ibnd,k,ispin)
                        pdos0_at(iatom)=pdos0_at(iatom)+bb
     &                        *weigd(ibnd,k,ispin)
                        if(ie1==ie2) then
                          pdos0_diag(ie1,li2,iatom)=
     &                      pdos0_diag(ie1,li2,iatom)+bb
     &                        *weigd(ibnd,k,ispin)
                        endif
                      enddo
                    enddo   !!! over ie1
                  enddo
                enddo   !!! over ie2
                l2m=l2m2
              enddo   !!! over li2
            enddo        !!! over iatom
            pdos0_at(0)=pdos0_at(0)+(1.d0-aa)*weigd(ibnd,k,ispin)
          enddo               !!! over bands
        enddo                  !!! over ind_k
        pdos_at(0,ispin)=pdos0_at(0)
        do isort=1,nsort
          pdos(:,isort,ispin)=0.d0
          pdos_at(isort,ispin)=0.d0
          pdos_diag(:,:,isort,ispin)=0.d0
          ndimr=nrel*lmb(isort)+1
          do iatom=1,natom
            if(is(iatom)==isort) then
              pdos_at(isort,ispin)=pdos_at(isort,ispin)+pdos0_at(iatom)
              do li2=1,ndimr
                pdos(li2,isort,ispin)=pdos(li2,isort,ispin)
     &	 	                       +pdos0(li2,iatom)
                pdos_diag(:,li2,isort,ispin)=
     &            pdos_diag(:,li2,isort,ispin)
     &	 	                       +pdos0_diag(:,li2,iatom)
              enddo
            endif
          enddo
          pdos_at(isort,ispin)=pdos_at(isort,ispin)/nhsort(isort)
          pdos(:,isort,ispin)=pdos(:,isort,ispin)/nhsort(isort)
          pdos_diag(:,:,isort,ispin)=pdos_diag(:,:,isort,ispin)
     &                              /nhsort(isort)
        enddo
      enddo   !!! over ispin
      deallocate(pdos0,pdos0_at,pdos0_diag)
      if(nproc_k/=1) then
        call DGOP(pdos,(nrel*maxb+1)*nsort*nspin,'  +',comm_k)
        call DGOP(pdos_at,(nsort+1)*nspin,'  +',comm_k)
        call DGOP(pdos_diag,maxntle*(nrel*maxb+1)*nsort*nspin,'  +',
     &            comm_k)
      endif
      END
