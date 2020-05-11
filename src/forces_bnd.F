      subroutine forces_bnd(f)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use units_mod
	  use solid_mod
      use vertex_mod
      implicit none
      real*8, intent(inout) :: f(3,nsort)
      integer :: isort,m,lm,ispin,n,l,m1,ii,iatom,ndimb,lm1,l1,lget,li1,
     &           i1,je,jj,km1,lme1,li,i,ie,km,lme,ind_k,k,nbas,j,j0,i0,
     &           ind,ndimb2,nbas0,ipiv(3),info
      real*8 :: pi2,pi2a,sb,kg(3),kg1(3),g1g(3),pi4om,kg22,v(3),
     &          fd,fermi_dirac,kg12,kg2,gg11,phase,de,aa(3,3)
      complex*16 :: cc
      integer, allocatable :: ind_y(:,:)
      real*8, allocatable :: f0(:,:),gb(:,:),yr(:,:),yi(:,:),a(:),b(:)
      complex*16, allocatable :: ff(:,:,:,:),oo(:,:,:,:),yc(:,:),
     &                           ts(:,:,:,:),os(:,:,:,:),tmp(:,:),
     &                           tmp1(:,:),tt(:,:,:),ss(:,:,:),
     &                           tmp2(:,:),rel(:,:),rel1(:,:)
      return
      pi2=pi+pi
      pi2a=pi2/par
      pi4om=(pi2+pi2)/amega
	  pv=0.d0
	  if(irel==2) pvj=0.d0
	  call vslli_mt(2,2,v_mt_h,v_mt_xc) ! H + XC
      allocate(f0(3,natom))
      f0=0.d0
      allocate(ff(2*limlb,2*limlb,natom,nspin))
      allocate(oo(2*limlb,2*limlb,natom,nspin))
      ff=(0.d0,0.d0)
      oo=(0.d0,0.d0)
      do ispin=1,nspin
        do iatom=1,natom
          isort=is(iatom)
          ndimb=nrel*(lmb(isort)+1)**2
          do lm1=1,ndimb
	        if(irel.le.1) then
	          l1=lget(lm1)
	          li1=l1+1
	        else
	          call getlimj(lm1,l1,i1,m1,li1,0)
	        endif
            je=n_e_apw(l1,isort)
            do jj=1,2
              km1=indbasa(jj,je,lm1,isort)
              if(jj==1) lme1=lm1
              if(jj==2) lme1=ndimb+lm1
              do lm=1,ndimb
	            if(irel.le.1) then
	              l=lget(lm)
	              li=l+1
	            else
	              call getlimj(lm,l,i,m,li,0)
	            endif
                ie=n_e_apw(l,isort)
                do ii=1,2
                  km=indbasa(ii,ie,lm,isort)
                  if(ii==1) lme=lm
                  if(ii==2) lme=ndimb+lm
                  if(lm==lm1) then
                    if(ii==1.and.jj==1) then
                      ff(lm,lm,iatom,ispin)=
     &                             ffhmt(1,1,ie,ie,li,isort,ispin)
                      oo(lm,lm,iatom,ispin)=
     &                             ffsmt(1,1,ie,ie,li,isort,ispin)
                    else if(ii==2.and.jj==1) then
                      ff(ndimb+lm,lm,iatom,ispin)=
     &                    0.5d0*(ffhmt(2,1,ie,ie,li,isort,ispin)
     &                          +ffhmt(1,2,ie,ie,li,isort,ispin))
                      oo(ndimb+lm,lm,iatom,ispin)=
     &                    ffsmt(2,1,ie,ie,li,isort,ispin)
                    else if(ii==1.and.jj==2) then
                      ff(lm,ndimb+lm,iatom,ispin)=
     &                    0.5d0*(ffhmt(2,1,ie,ie,li,isort,ispin)
     &                          +ffhmt(1,2,ie,ie,li,isort,ispin))
                      oo(lm,ndimb+lm,iatom,ispin)=
     &                    ffsmt(1,2,ie,ie,li,isort,ispin)
                    else if(ii==2.and.jj==2) then
                      ff(ndimb+lm,ndimb+lm,iatom,ispin)=
     &                             ffhmt(2,2,ie,ie,li,isort,ispin)
                      oo(ndimb+lm,ndimb+lm,iatom,ispin)=
     &                             ffsmt(2,2,ie,ie,li,isort,ispin)
                    endif
                  endif
                  ff(lme,lme1,iatom,ispin)=ff(lme,lme1,iatom,ispin)
     &                                    +pv(km,km1,iatom,ispin)
                  if(irel==2) ff(lme,lme1,iatom,ispin)=
     &                          ff(lme,lme1,iatom,ispin)
     &                          +(0.d0,1.d0)*pvj(km,km1,iatom)
                enddo  !! over ii
              enddo  !! over lm
            enddo  !! over jj
          enddo   !! over lm1
        enddo   !! over iatom
      enddo !! over ispin
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        nbas0=nbask(k)
        nbas=nbas0/nrel
        allocate(gb(3,nbas))
        do i=1,nbas
          gb(:,i)=gbs(:,indgb(i,k))
        enddo
        allocate(yc(max_bs,nbas0))
        allocate(yr(max_bs,nbas0))
        allocate(yi(max_bs,nbas0))
        yr=0.d0
        yi=0.d0
        allocate(ind_y(2*limlb,natom))
        call setk2(pnt(1,k),gb,yr,yi,nbas0,ind_y)
        yc=dcmplx(yr,yi)
        deallocate(yr,yi)
        allocate(ts(nbas0,nbas0,3,natom))
        allocate(os(nbas0,nbas0,3,natom))
        ts=(0.d0,0.d0)
        os=(0.d0,0.d0)
        if(irel==2) then
          allocate(rel(nbas0,nbas0))
          allocate(rel1(nbas0,nbas0))
          call forces_rel_bnd(rel,rel1,gb,pnt(1,k),nbas0,indgb(1,k))
        endif
        deallocate(gb)
        do j=1,nbas0
          if(j<=nbas) jj=j
          if(j>nbas) jj=j-nbas
          j0=indgb(jj,k)
          kg1=pi2a*(pnt(:,k)+gbs(:,j0))
          kg12=dot_product(kg1,kg1)
          do i=1,nbas0
            if(i<=nbas) ii=i
            if(i>nbas) ii=i-nbas
            i0=indgb(ii,k)
            if(i0==j0) cycle
            kg=pi2a*(pnt(:,k)+gbs(:,i0))
            kg2=dot_product(kg,kg)
            kg22=0.5d0*(kg12+kg2)
            g1g=kg1-kg
            gg11=sqrt(dot_product(g1g,g1g))
            v=g1g/gg11
            do iatom=1,natom
              isort=is(iatom)
              sb=gg11*smt(isort)
              sb=(sin(sb)/sb-cos(sb))/sb
              sb=smt(isort)**2*sb
              phase=par*dot_product(g1g,tau(:,iatom))
              cc=(0.d0,1.d0)*dcmplx(cos(phase),sin(phase))*sb
              if(irel<2) then
                os(i,j,:,iatom)=pi4om*cc*v
                ts(i,j,:,iatom)=kg22*os(i,j,:,iatom)
              else
                os(i,j,:,iatom)=pi4om*cc*v*rel(i,j)
                ts(i,j,:,iatom)=pi4om*cc*v*rel1(i,j)
              endif
            enddo
          enddo
        enddo
        allocate(tt(nbas0,nbas0,3))
        allocate(ss(nbas0,nbas0,3))
        do ispin=1,nspin
          n=n_bnd(k,ispin)
          allocate(a(n))
          allocate(b(n))
          do iatom=1,natom
            isort=is(iatom)
            ndimb=nrel*(lmb(isort)+1)**2
            ind=ind_y(1,iatom)
            ndimb2=2*ndimb
            tt=ts(:,:,:,iatom)
            ss=os(:,:,:,iatom)
            allocate(tmp(nbas0,ndimb2))
            allocate(tmp1(nbas0,nbas0))
            allocate(tmp2(nbas0,nbas0))
            call zgemm('c','n',nbas0,ndimb2,ndimb2,(1.d0,0.d0),
     &                 yc(ind,1),max_bs,ff(1,1,iatom,ispin),2*limlb,
     &                 (0.d0,0.d0),tmp,nbas0)
            call zgemm('n','n',nbas0,nbas0,ndimb2,(1.d0,0.d0),tmp,nbas0,
     &                 yc(ind,1),max_bs,(0.d0,0.d0),tmp1,nbas0)
            call zgemm('c','n',nbas0,ndimb2,ndimb2,(1.d0,0.d0),
     &                 yc(ind,1),max_bs,oo(1,1,iatom,ispin),2*limlb,
     &                 (0.d0,0.d0),tmp,nbas0)
            call zgemm('n','n',nbas0,nbas0,ndimb2,(1.d0,0.d0),tmp,nbas0,
     &                 yc(ind,1),max_bs,(0.d0,0.d0),tmp2,nbas0)
            deallocate(tmp)
            do j=1,nbas0
              if(j<=nbas) jj=j
              if(j>nbas) jj=j-nbas
              j0=indgb(jj,k)
              do i=1,nbas0
                if(i<=nbas) ii=i
                if(i>nbas) ii=i-nbas
                i0=indgb(ii,k)
                if(i0==j0) cycle
                g1g=pi2a*(gbs(:,j0)-gbs(:,i0))
                cc=(0.d0,1.d0)*tmp1(i,j)
                tt(i,j,:)=tt(i,j,:)-g1g*cc
                cc=(0.d0,1.d0)*tmp2(i,j)
                ss(i,j,:)=ss(i,j,:)-g1g*cc
              enddo
            enddo
            deallocate(tmp1,tmp2)
            allocate(tmp(nbas0,n))
            do i=1,3
              call zgemm('n','n',nbas0,n,nbas0,(1.d0,0.d0),tt(1,1,i),
     &                   nbas0,ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &                   (0.d0,0.d0),tmp,nbas0)
              do j=1,n
                a(j)=dot_product(ev_bnd(1:nbas0,j,ind_k,ispin),tmp(:,j))
              enddo
              call zgemm('n','n',nbas0,n,nbas0,(1.d0,0.d0),ss(1,1,i),
     &                   nbas0,ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &                   (0.d0,0.d0),tmp,nbas0)
              do j=1,n
                b(j)=dot_product(ev_bnd(1:nbas0,j,ind_k,ispin),tmp(:,j))
              enddo
              do j=1,n
                de=e_bnd(j,k,ispin)-chem_pot
                fd=fermi_dirac(de)*wgt(k)
                f0(i,iatom)=f0(i,iatom)+fd*(a(j)-e_bnd(j,k,ispin)*b(j))
              enddo
            enddo
            deallocate(tmp)
          enddo   !!! over iatom
          deallocate(a,b)
        enddo   !! over ispin
        deallocate(yc,ts,os,tt,ss,ind_y)
        if(irel==2) deallocate(rel,rel1)
      enddo   !! over ind_k
      deallocate(ff,oo)
      if(nspin==1.and.irel/=2) f0=f0*2
	if(nproc_k/=1) call dgop(f0,3*natom,'  +',comm_k)
      call sym_force(f0)
      do isort=1,nsort
        iatom=iat_1(isort)
        f(:,isort)=f(:,isort)+f0(:,iatom)
      enddo
      deallocate(f0)
      if(maswrk) then
        write(iun,*)' Forces after BND'
        do isort=1,nsort
          write(iun,'(a4,3(1x,f12.5))')txtel(isort),(f(i,isort),i=1,3)
        enddo
        write(iun,*)' Forces in Fractional Coordinates'
        do isort=1,nsort
          v=f(:,isort)
          aa=rbas
          call dgesv(3,1,aa,3,ipiv,v,3,info)
          write(iun,'(a4,3(1x,f12.5))')txtel(isort),(v(i),i=1,3)
        enddo
      endif
      end
