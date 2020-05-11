      subroutine orbitals_red
	  use atom_mod
	  use heg_mod
	  use manager_mod
	  use parallel_mod
      use solid_mod
	  use units_mod
	  use vertex_mod
	  implicit none
	  integer :: i,isort,l,nlm,li,ie,ii,jj,lm,lf,lfm,iend,mj,iatom
      if(.not.allocated(lfun_red)) allocate(lfun_red(nsort))
      if(.not.allocated(lfunm_red)) allocate(lfunm_red(nsort))
	  maxel_red=0
	  maxlfun_red=0
	  nfun_red=0
      do isort=1,nsort
	    lfun_red(isort)=0
	    lfunm_red(isort)=0
	    do li=1,nrel*lmb_red(isort)+1
          if(irel.ne.2) then
		    l=li-1
            nlm=l+l+1
          else if(irel.eq.2) then
		    l=li/2
            nlm=((li+1)/2)*2
	      endif
	      do ie=1,ntle(l,isort)
            if(correlated(ie,l,isort)=='N') cycle
	        lfun_red(isort)=lfun_red(isort)+1
	        lfunm_red(isort)=lfunm_red(isort)+nlm
            if(correlated(ie,l,isort)=='A') cycle
	        if(augm(ie,l,isort)=='LAW') then
	          lfun_red(isort)=lfun_red(isort)+1
	          lfunm_red(isort)=lfunm_red(isort)+nlm
	        endif
            if(augm(ie,l,isort)=='AWL') then
              lfun_red(isort)=lfun_red(isort)+1
              lfunm_red(isort)=lfunm_red(isort)+nlm
            endif
	      enddo
	    enddo
	    maxel_red=max(maxel_red,lfunm_red(isort))
	    maxlfun_red=max(maxlfun_red,lfun_red(isort))
	    nfun_red=nfun_red+lfunm_red(isort)*nhsort(isort)
      enddo
c -------------------------------------------------------------------
      if(.not.allocated(le_red)) allocate(le_red(maxlfun_red,nsort))
      if(.not.allocated(lme_red)) allocate(lme_red(maxel_red,nsort))
      if(.not.allocated(lf_isz_red)) 
     &  allocate(lf_isz_red(maxel_red,nsort))
      if(.not.allocated(indbasa_red)) 
     &  allocate(indbasa_red(2,maxntle,limlb_red,nsort))
      indbasa_red=0
      lf_isz_red=0
      do isort=1,nsort        
	    lf=0
	    lfm=0
        do li=1,nrel*lmb_red(isort)+1
	      if(irel.lt.2) then
	        l=li-1
	        jj=l
	      else
	        l=li/2
	        ii=mod(li,2)
	        if(ii.eq.0) ii=-1
            jj=l+l+ii
	      endif
          do ie=1,ntle(l,isort)
            if(correlated(ie,l,isort)=='N') cycle
            iend=1
            if(augm(ie,l,isort)=='LAW') iend=2
            if(augm(ie,l,isort)=='AWL') iend=2
            if(correlated(ie,l,isort)=='A') iend=1
            do i=1,iend
	          lf=lf+1
	          le_red(lf,isort)=indfun0(i,ie,li,isort)
              do mj=-jj,jj,nrel
                lfm=lfm+1
                lf_isz_red(lfm,isort)=lf
                if(irel.lt.2) then
			      lm=l*(l+1)+mj+1
	            else
		          call getlimj(lm,l,ii,mj,li,1)
	            endif
	            indbasa_red(i,ie,lm,isort)=lfm
	            lme_red(lfm,isort)=indbasa(i,ie,lm,isort)
	          enddo
	        enddo
	      enddo
	    enddo
	  enddo
c ---------------------------------------------------------
      if(.not.allocated(io_lem_red)) allocate(io_lem_red(natom))
      if(.not.allocated(io_lem2_red)) allocate(io_lem2_red(natom))
	  jj=1
	  ii=1
	  do iatom=1,natom
	    isort=is(iatom)
	    io_lem_red(iatom)=jj
	    io_lem2_red(iatom)=ii
	    jj=jj+lfunm_red(isort)
	    ii=ii+lfunm_red(isort)**2
	  enddo
	  n_lem2_red=ii-1
      end
