      subroutine seteny
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
      logical :: low
      integer :: ispin,isort,ndimr,li,l,ie,lm,ndim,m,k,ibnd,iatom,
     &	         km,n_bnd0,n_bnd1,ind_k,in,jn,lm0,n,i,ii
      real*8 :: tos0,tos1,admixe,enynew,tm,adm
      integer, allocatable :: ibm(:,:)
      real*8, allocatable :: eb(:,:),t0(:),t1(:)
	  admixe=0.5d0
	  adm=0.5d0
	  if(n_semic/=0) then
	    allocate(eb(2,n_semic))
	    allocate(ibm(2,n_semic))
	  endif
      do ispin=1,nspin
        if(n_semic/=0) then
          eb(1,:)=1.d10
          eb(2,:)=-1.d10
          do k=1,npnt
            do ie=1,n_semic
              eb(1,ie)=min(eb(1,ie),e_bnd(ie,k,ispin))
              eb(2,ie)=max(eb(2,ie),e_bnd(ie,k,ispin))
            enddo
          enddo
          n=1
          ibm(1,n)=1
          do ie=1,n_semic-1
            if(eb(2,ie)<eb(1,ie+1)) then
              ibm(2,n)=ie
              n=n+1
              ibm(1,n)=ie+1
            endif
          enddo
          ibm(2,n)=n_semic
          allocate(t0(n))
          allocate(t1(n))
        endif
        do isort=1,nsort
          ndimr=nrel*lmb(isort)+1
          lm0=0
	      do li=1,ndimr
	        if(irel.le.1) l=li-1
	        if(irel.eq.2) l=li/2
            if(irel.ne.2) ndim=li+li-1
            if(irel.eq.2) ndim=((li+1)/2)*2
	        do ie=1,ntle(l,isort)
			  in=1
			  if(augm(ie,l,isort)=='LAW') in=2
              if(augm(ie,l,isort)=='AWL') in=2
			  low=.false.
              if(augm(ie,l,isort)=='LOC'.and.ie<n_e_apw(l,isort))
     &          low=.true.
              if(low) then
                t0=0.d0
                t1=0.d0
              endif
	          tos0=0.d0
	          tos1=0.d0
	          do ind_k=1,ndim3_k(me_k+1)
	            k=n3_mpi_k(me_k+1)+ind_k
	            if(low) then
			      n_bnd0=1
			      n_bnd1=n_semic
	            else if(idmd(ie,l,isort)/=5) then
			      n_bnd0=n_semic+1
			      n_bnd1=n_bnd(k,ispin)
                else if(idmd(ie,l,isort)==5) then
                  n_bnd0=n_semic+1
                  n_bnd1=n_semic+nbndf_bnd_0(2)
	            endif
	            do ibnd=n_bnd0,n_bnd1
                  if(low) then
                    ii=0
                    do i=1,n
                      if(ibm(1,i)<=ibnd.and.ibm(2,i)>=ibnd) ii=i
                    enddo
                  endif
	              tm=0.d0
	              do iatom=1,natom
	                if(is(iatom)==isort) then
	                  do jn=1,in
	                    lm=lm0
	                    do m=1,ndim
	                      lm=lm+1
	                      km=io_lem(iatom)+indbasa(jn,ie,lm,isort)-1
	                      tm=tm+conjg(z_bnd(km,ibnd,ind_k,ispin))
     &		 		                 *z_bnd(km,ibnd,ind_k,ispin)
                        enddo
	                  enddo
	                endif
	              enddo
	              tm=tm/nhsort(isort)*wgt(k)
                  if(low) then
                    t0(ii)=t0(ii)+g_full_00(ibnd,ind_k,ispin)*tm
                    t1(ii)=t1(ii)+e_bnd(ibnd,k,ispin)
     &                       *g_full_00(ibnd,ind_k,ispin)*tm
                  else if(idmd(ie,l,isort)/=5) then
 	                tos0=tos0+g_full_00(ibnd,ind_k,ispin)*tm
	                tos1=tos1+e_bnd(ibnd,k,ispin)
     &			             *g_full_00(ibnd,ind_k,ispin)*tm
                  else if(idmd(ie,l,isort)==5) then
                    tos0=tos0+tm
                    tos1=tos1+e_bnd(ibnd,k,ispin)*tm
                  endif
	            enddo   !! over ibnd
	          enddo   !! over ind_k
	          if(nproc_k/=1) then
                if(low) then
                  call DGOP(t0,n,'  +',comm_k)
                  call DGOP(t1,n,'  +',comm_k)
                else
                  call DGOP(tos0,1,'  +',comm_k)
                  call DGOP(tos1,1,'  +',comm_k)
                endif
	          endif
              if(idmd(ie,l,isort)==0.or.idmd(ie,l,isort)==5) then
                if(low) then
                  tos0=0.d0
                  do i=1,n
                    if(t0(i)>tos0) then
                      tos0=t0(i)
                      ii=i
                    endif
                  enddo
                  tos1=t1(ii)
                endif
	            if(abs(tos0).gt.1.d-8) then
                  enynew=tos1/tos0
                  enynew=admixe*enynew
     &			        +(1.d0-admixe)*eny(ie,li,isort,ispin)
                  eny(ie,li,isort,ispin)=enynew
	            endif
	          endif
              if(idmd(ie,l,isort)==4) then
                enynew=chem_pot
                enynew=adm*enynew+(1.d0-adm)*eny(ie,li,isort,ispin)
                eny(ie,li,isort,ispin)=enynew
	          endif
	        enddo   !! over ie
	        lm0=lm0+ndim
	      enddo
	    enddo
	    if(n_semic/=0) deallocate(t0,t1)
	  enddo
	  if(n_semic/=0) deallocate(eb,ibm)
      end
