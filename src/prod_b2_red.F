      subroutine prod_b2_red(ndim0,key,nl,ind_l,limax)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
c     key=0 - Run to find the dimensions
c     key=1 - Run to fill the arrays
      integer, intent(in) :: ndim0,limax,nl(limax,limax),
     &                       ind_l(limax,limax,maxpb+1),key
      integer :: ii,isort,n1,n2,l2,n0,ispin,in,li,mt,ir,l,jn,li1,mt1,l1,
     &	       nl0,j,i,lwork,il,iu,nev,info,n,i0,m2,lm2,ii2,iatom,
     &           in0,jn0,l2max,ie,je,jspin
      real*8 :: c2,s00,dqdall,abstol,vl,vu,sqe,s,ac
      integer, allocatable :: iwork(:),ifail(:)
      real*8, allocatable :: s0(:,:),evec(:,:),eval(:),pp(:,:),pp1(:,:),
     &	                   qq(:,:),qq1(:,:),pb(:,:),work(:),wrk(:),
     &                       dpb(:,:)
      allocate(s0(ndim0,ndim0),evec(ndim0,ndim0),eval(ndim0),
     &	     pp(0:maxnrad,2),pp1(0:maxnrad,2),qq(0:maxnrad,2),
     &         qq1(0:maxnrad,2),pb(0:maxnrad,ndim0),work(0:maxnrad))
      allocate(dpb(0:maxnrad,ndim0))
      c2=clight**2
      if(key.eq.0) then
        n_pbmt_red=0
        maxpbr_red=0
        max_e_pb_red=0
      endif
      if(key==1) then
        fun_pb_red=0.d0
        dfun_pb_red=0.d0
      endif
      l2max=maxval(ind_l)
      ii=1
      do isort=1,nsort
        if(key==0.and.maswrk) write(iun,*)
     &	  ' Reduced Product Basis for isort ',isort
        if(key.eq.0) then
          lfun_pb_red(isort)=0
          n_pbmt0_red(isort)=0
        endif
        if(key.eq.1) then
          n1=0
          n2=0
        endif
        do l2=0,min(l2max,lmpb_red(isort))
          if(key==0.and.maswrk) write(iun,*)' Eigen values for L = ',l2
            n0=0
            do ispin=1,nspin
              do in0=1,lfun_red(isort)
                in=le_red(in0,isort)
                li=li_val(in,isort)
                if(irel.le.1) l=li-1
                if(irel.eq.2) l=li/2
                ie=ie_val(in,isort)
                if(correlated(ie,l,isort)=='A') then
	              do ir=1,nrad(isort)
	                pp(ir,1)=p_f(ir,ie,li,isort,ispin)/r(ir,isort)
	                pp(ir,2)=p_f(ir,ie,li,isort,ispin)/r(ir,isort)
	                if(irel.ge.1) qq(ir,1)=q_f(ir,ie,li,isort,ispin)
     &                                    /r(ir,isort)
	                if(irel.ge.1) qq(ir,2)=q_f(ir,ie,li,isort,ispin)
     &                                    /r(ir,isort)
                  enddo
                  call fit_zero(pp(0,1),r(0,isort))
                  if(irel>0) call fit_zero(qq(0,1),r(0,isort))
                else if(correlated(ie,l,isort)=='L') then
	              mt=ind_wf(in,isort)
	              do ir=0,nrad(isort)
	                pp(ir,1)=gfun(mt+ir,ispin)
	                if(irel.ge.1) qq(ir,1)=gfund(mt+ir,ispin)
	                pp(ir,2)=gfun1(mt+ir,ispin)
	                if(irel.ge.1) qq(ir,2)=gfund1(mt+ir,ispin)
                  enddo
                endif
c	          do jspin=1,ispin
                jspin=ispin
	            do jn0=1,lfun_red(isort)
	              jn=le_red(jn0,isort)
			      li1=li_val(jn,isort)
	              if(irel.le.1) l1=li1-1
	              if(irel.eq.2) l1=li1/2
                  je=ie_val(jn,isort)
                  if(correlated(je,l1,isort)=='A') then
	                do ir=1,nrad(isort)
	                  pp1(ir,1)=p_f(ir,je,li1,isort,jspin)/r(ir,isort)
                      if(irel.ge.1) qq1(ir,1)=q_f(ir,je,li1,isort,jspin)
     &                                       /r(ir,isort)
	                  pp1(ir,2)=p_f(ir,je,li1,isort,jspin)/r(ir,isort)
                      if(irel.ge.1) qq1(ir,2)=q_f(ir,je,li1,isort,jspin)
     &                                       /r(ir,isort)
                    enddo
                    call fit_zero(pp1(0,1),r(0,isort))
                    if(irel>0) call fit_zero(qq1(0,1),r(0,isort))
                  else if(correlated(je,l1,isort)=='L') then
	                mt1=ind_wf(jn,isort)
	                do ir=0,nrad(isort)
	                  pp1(ir,1)=gfun(mt1+ir,jspin)
	                  if(irel.ge.1) qq1(ir,1)=gfund(mt1+ir,jspin)
	                  pp1(ir,2)=gfun1(mt1+ir,jspin)
	                  if(irel.ge.1) qq1(ir,2)=gfund1(mt1+ir,jspin)
                    enddo
                  endif
	              do nl0=1,nl(li,li1)
	                if(ind_l(li,li1,nl0).eq.l2) then
	                  n0=n0+1
	                  do ir=0,nrad(isort)
	                    pb(ir,n0)=pp(ir,1)*pp1(ir,1)
                        dpb(ir,n0)=pp(ir,2)*pp1(ir,1)+pp(ir,1)*pp1(ir,2)
                      enddo
                      if(irel/=0) then
	                    do ir=0,nrad(isort)
	                      pb(ir,n0)=pb(ir,n0)+qq(ir,1)*qq1(ir,1)/c2
	                      dpb(ir,n0)=dpb(ir,n0)+(qq(ir,2)*qq1(ir,1)
     &                              +qq(ir,1)*qq1(ir,2))/c2
                        enddo
                      endif
c ---- Make it orthogonal to the constant ------------------------
                      if(mode_pb==1.and.l2==0) then
                        do ir=0,nrad(isort)
                          work(ir)=pb(ir,n0)*dr(ir,isort)
     &                                      *r(ir,isort)**2
                        enddo
                        s00=dqdall(h(isort),work,nrad(isort))
                        do ir=0,nrad(isort)
                          pb(ir,n0)=pb(ir,n0)-ac*ac*s00
                        enddo
                      endif
                      do ir=0,nrad(isort)
                        work(ir)=pb(ir,n0)*pb(ir,n0)*dr(ir,isort)
     &					  	                        *r(ir,isort)**2
                      enddo
                      s00=dqdall(h(isort),work,nrad(isort))
	                  s00=1.d0/sqrt(s00)
	                  call dscal(nrad(isort)+1,s00,pb(0,n0),1)
	                  call dscal(nrad(isort)+1,s00,dpb(0,n0),1)
	                endif
	              enddo
                enddo   !! over jn
c              enddo  !! over jspin
	        enddo   !! over in
	      enddo   !! over ispin
	      do j=1,n0
	        do i=1,j
              do ir=0,nrad(isort)
                pp(ir,1)=pb(ir,i)*pb(ir,j)*dr(ir,isort)*r(ir,isort)**2
              enddo
              s0(i,j)=dqdall(h(isort),pp,nrad(isort))
	        enddo
	      enddo
c ---------- Diagonalize ----------------------------------------
	      abstol=1.d-70
          lwork=max(1,8*n0)
          allocate(wrk(lwork),iwork(5*n0),ifail(n0))
	      call DSYEVX('V','A','U',n0,s0,ndim0,VL,VU,IL,IU,ABSTOL,nev,
     &		          eval,evec,ndim0,wrk,lwork,iwork,ifail,INFO) 
	      if(info.ne.0) then
	        if(maswrk) write(iun,*)'Prod_b2: Info = ',info
	        if(info.gt.n0) call ending
	      endif
          deallocate(wrk,iwork,ifail)
	      n=0
          if(mode_pb==1.and.l2==0) n=1
          if(lim_pb_mt_red(l2,isort)/=0) then
	        do i=n0,1,-1
	          if(eval(i)/eval(n0).ge.eps_pb_vrt) n=n+1
	          if(n==lim_pb_mt_red(l2,isort)) exit
            enddo
          endif
	      if(n/=0.and.key.eq.0) then
	        ntle_pb_red(l2,isort)=n
	        lfun_pb_red(isort)=lfun_pb_red(isort)+n
	        n_pbmt0_red(isort)=n_pbmt0_red(isort)+n*(l2+l2+1)
	        max_e_pb_red=max(n,max_e_pb_red)
	        if(maswrk) then
              if(mode_pb==1.and.l2==0) then
                write(iun,'(7(1x,f9.6))')1.d0,(eval(i)/eval(n0),
     &                                   i=n0,n0-n+2,-1)
              else
                write(iun,'(7(1x,f9.6))')(eval(i)/eval(n0),
     &                                   i=n0,n0-n+1,-1)
              endif
	        endif
	      else if(key.eq.1) then
	        do i=1,n
	          n1=n1+1
	          ind_pb_red(n1,isort)=ii
	          lval_pb_red(n1,isort)=l2
              ind_prod0_red(i,l2,isort)=n1
              if(mode_pb==1.and.l2==0.and.i==1) then
                do ir=0,nrad(isort)
                  fun_pb(ii+ir)=ac
                  dfun_pb(ii+ir)=0.d0
                enddo
              else
                i0=n0-i+1
                if(mode_pb==1.and.l2==0) i0=i0+1
                sqe=sqrt(eval(i0))
                do j=1,n0
                  s=evec(j,i0)/sqe
                  do ir=0,nrad(isort)
                    fun_pb_red(ii+ir)=fun_pb_red(ii+ir)+pb(ir,j)*s
                    dfun_pb_red(ii+ir)=dfun_pb_red(ii+ir)+dpb(ir,j)*s
                  enddo
                enddo
              endif
              ii=ii+nrad(isort)+1
              do m2=-l2,l2
                n2=n2+1
                lm2=l2*(l2+1)+m2+1
                lm_pbmt_red(n2,isort)=lm2
                lf_pb_red(n2,isort)=n1
                ind_prod_red(i,lm2,isort)=n2
              enddo
            enddo
          endif
        enddo   !! over l2
        if(key.eq.0) then
          maxlfpb_red=max(maxlfpb_red,lfun_pb_red(isort))
          n_pbmtm_red=max(n_pbmtm_red,n_pbmt0_red(isort))
          n_pbmt_red=n_pbmt_red+n_pbmt0_red(isort)*nhsort(isort)
          maxpbr_red=maxpbr_red+lfun_pb_red(isort)*(nrad(isort)+1)
        endif
      enddo   !! over isort
      if(key.eq.0) then
        ii2=1
        do iatom=1,natom
          isort=is(iatom)
          iopb_red(iatom)=ii2
          ii2=ii2+n_pbmt0_red(isort)
        enddo
        ii2=1
        do isort=1,nsort
          iopb0_red(isort)=ii2
          ii2=ii2+n_pbmt0_red(isort)
        enddo
        if(maswrk) then
          do isort=1,nsort
            write(iun,*)' Reduced Product Basis info for ',txtel(isort)
            write(iun,'(1x,a1,4x,a4,1x,a11)')'L','ntle','ntle*(2L+1)'
            do l2=0,min(l2max,lmpb_red(isort))
              write(iun,'(i2,5x,i3,8x,i4)')l2,ntle_pb_red(l2,isort),
     &		                             (l2+l2+1)*ntle_pb_red(l2,isort)
            enddo
          enddo
        endif
      endif
      deallocate(s0,evec,eval,pp,pp1,qq,qq1,pb,work,dpb)
      if(key==1) call timel('**** prod_b2_red finished **********')
      end
