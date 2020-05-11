      subroutine bnd_curv
	use atom_mod
	use manager_mod
	use parallel_mod
	use units_mod
	use solid_mod
      implicit none
	integer :: k,ispin,n,alf,bet,ib,ie,i,j,ind,ibl,ibnd,i01,nbl,i1,i0,
     &           ndeg,j0,ind0
      real*8 :: dd,q(3),q2,q1,qsm,e1d(3),e2d(3),cv(3,3)
      complex*16 :: cc,c2(3,3),e1(3,3),e2(3,3),tmp(3,3),vv(3)
      real*8, allocatable :: grd(:,:),crv(:,:,:),ce(:,:)
      complex*16, allocatable :: ea(:,:,:),eab(:,:,:),c0(:)
      return
      k=1
      ispin=1
      q(1)=-1.d0
      q(2)=1.d0
      q(3)=1.d0
      qsm=0.001d0
      q2=dot_product(q,q)
      q1=qsm/sqrt(q2)
      q=q*q1
      q2=q2*q1**2
      n=n_bnd(k,ispin)
      allocate(ea(n,n,3))
      call bands_gradient(nbndf,n,k,ispin,z_bnd(1,1,k,ispin),
     &                    ev_bnd(1,1,k,ispin),ea)
      ea=(0.d0,-2.d0)*ea
      allocate(eab(k_irrep_max,3,3))
      do ibl=1,n_block_k(k,ispin)
        ie=last_block_k(ibl,k,ispin)
        if(ibl==1) ib=1
        if(ibl/=1) ib=last_block_k(ibl-1,k,ispin)+1
        ind=entry_block_k(ibl,k,ispin)-1
        do j=ib,ie
          do i=ib,ie
            ind=ind+1
            do bet=1,3
              do alf=1,3
                cc=(0.d0,0.d0)
                if(alf==bet.and.i==j) cc=(1.d0,0.d0)
                do ibnd=1,n
                  if(ibnd>=ib.and.ibnd<=ie) cycle
                  dd=e_bnd(ib,k,ispin)-e_bnd(ibnd,k,ispin)
                  cc=cc+ea(i,ibnd,alf)*ea(ibnd,j,bet)/dd
                enddo
                eab(ind,alf,bet)=2.d0*cc
              enddo
            enddo
          enddo
        enddo
      enddo
c ------------------- Form C0 -------------------------------
      allocate(c0(k_irrep_max))
      c0=(0.d0,0.d0)
      do ibl=1,n_block_k(k,ispin)
        ie=last_block_k(ibl,k,ispin)
        if(ibl==1) ib=1
        if(ibl/=1) ib=last_block_k(ibl-1,k,ispin)+1
        nbl=ie-ib+1
        ind=entry_block_k(ibl,k,ispin)
        if(nbl==1) then
          c0(ind)=(1.d0,0.d0)
        else
          e1=(0.d0,0.d0)
          dd=0.d0
          do j=ib,ie
            do i=ib,ie
              do alf=1,3
                e1(i-ib+1,j-ib+1)=e1(i-ib+1,j-ib+1)+q(alf)*ea(i,j,alf)
              enddo
              dd=max(dd,abs(e1(i-ib+1,j-ib+1)))
            enddo
          enddo
          e1d=0.d0
          if(dd>qsm*1.d-8) call eig_val_solver(nbl,3,e1,e1d)
          ind=ind-1
          do j=ib,ie
            do i=ib,ie
              ind=ind+1
              if(dd>qsm*1.d-8) c0(ind)=e1(i-ib+1,j-ib+1)
              if(dd<=qsm*1.d-8.and.j==i) c0(ind)=(1.d0,0.d0)
            enddo
          enddo
c ------------------------ E2 ------------------------------------
          e2=(0.d0,0.d0)
          ind=entry_block_k(ibl,k,ispin)-1
          do j=ib,ie
            do i=ib,ie
              ind=ind+1
              do alf=1,3
                do bet=1,3
                  e2(i-ib+1,j-ib+1)=e2(i-ib+1,j-ib+1)+q(alf)*q(bet)
     &                             *eab(ind,alf,bet)
                enddo
              enddo
            enddo
          enddo
          e2=e2/2
          ind=entry_block_k(ibl,k,ispin)
          call zgemm('n','n',nbl,nbl,nbl,(1.d0,0.d0),e2,3,c0(ind),nbl,
     &               (0.d0,0.d0),tmp,3)
          call zgemm('c','n',nbl,nbl,nbl,(1.d0,0.d0),c0(ind),nbl,tmp,3,
     &               (0.d0,0.d0),e2,3)
c -------------------------------------------------------------   
	    i1=0
1	    i0=i1+1
	    ndeg=1
          i01=i0+1
	    do i=i01,nbl
	      if(abs(e1d(i)-e1d(i0))<1.d-8) ndeg=ndeg+1
          enddo
          if(ndeg/=1) then
            call eig_val_solver(ndeg,3,e2(i0,i0),e2d)
            if(ibl==4) print *,e2d(1),e2d(2),e2d(3)
            ind=entry_block_k(ibl,k,ispin)-1
            do j=1,nbl
              do i=1,nbl
                ind=ind+1
                c2(i,j)=c0(ind)
              enddo
            enddo
            ind=entry_block_k(ibl,k,ispin)+nbl*(i0-1)
            call zgemm('n','n',nbl,ndeg,ndeg,(1.d0,0.d0),c2(1,i0),3,
     &                 e2(i0,i0),3,(0.d0,0.d0),c0(ind),3)
          endif
	    i1=i0+ndeg-1
          if(i1<nbl) goto 1
        endif
      enddo
      allocate(grd(3,n))
      allocate(crv(3,3,n))
      grd=0.d0
      crv=0.d0
      allocate(ce(3,n))
      do ibnd=1,n
        ibl=ind_block_k(ibnd,k,ispin)
        ie=last_block_k(ibl,k,ispin)
        if(ibl==1) ib=1
        if(ibl/=1) ib=last_block_k(ibl-1,k,ispin)+1
        nbl=ie-ib+1
        ind=entry_block_k(ibl,k,ispin)+nbl*(ibnd-ib)
        vv(1:nbl)=c0(ind:ind+nbl-1)
        ind=entry_block_k(ibl,k,ispin)
        do j=ib,ie
          j0=j-ib+1
          do i=ib,ie
            i0=i-ib+1
            ind0=ind+nbl*(j0-1)+i0-1
            do alf=1,3
              grd(alf,ibnd)=grd(alf,ibnd)+conjg(vv(i0))
     &                                   *ea(i,j,alf)*vv(j0)
              do bet=1,3
                cc=(eab(ind0,alf,bet)+eab(ind0,bet,alf))/2
                crv(alf,bet,ibnd)=crv(alf,bet,ibnd)+conjg(vv(i0))
     &                                             *cc*vv(j0)
              enddo
            enddo
          enddo
        enddo
        cv=crv(:,:,ibnd)
        call eig_val_solver_r(3,3,cv,ce(1,ibnd))
      enddo
      deallocate(ea,eab,c0,grd,crv,ce)
      end