      subroutine rpb_from_fpb_1(isort,jsort,nr0,nf0,pr,pf)
	use atom_mod
	use vertex_mod
	implicit none
	integer, intent(in) :: isort,jsort,nr0,nf0
	real*8, intent(in) :: pf(nf0,nf0)
	real*8, intent(out) :: pr(nr0,nr0)
	integer :: lget,nlm,lm1,lm,l1,l,ie,ie1,ie2,ie3,km,km1,km2,km3,nlm1
	nlm=(lmpb_red(isort)+1)**2
	nlm1=(lmpb_red(jsort)+1)**2
	pr=0.d0
	do lm1=1,nlm1
	  l1=lget(lm1)
	  do ie3=1,ntle_pb(l1,jsort)
	    km3=ind_prod(ie3,lm1,jsort)
	    do ie1=1,ntle_pb_red(l1,jsort)
	      km1=ind_prod_red(ie1,lm1,jsort)
	      do lm=1,nlm
	        l=lget(lm)
	        do ie2=1,ntle_pb(l,isort)
	          km2=ind_prod(ie2,lm,isort)
	          do ie=1,ntle_pb_red(l,isort)
	            km=ind_prod_red(ie,lm,isort)
	            pr(km,km1)=pr(km,km1)+rpb_fpb(ie,ie2,l,isort)
     &				      *pf(km2,km3)*rpb_fpb(ie1,ie3,l1,jsort)
	          enddo
	        enddo
	      enddo
	    enddo
	  enddo
	enddo
      end