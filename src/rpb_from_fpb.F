      subroutine rpb_from_fpb(pl,pf)
	use atom_mod
	use solid_mod
	use vertex_mod
	implicit none
	real*8, intent(in) :: pf(n_pbmt,n_pbmt)
	real*8, intent(out) :: pl(n_pbmt_red,n_pbmt_red)
	integer :: lget,ie,ie1,ie2,ie3,km,km1,km2,km3,iatom,jatom,isort,
     &	         jsort,j_f,j_l,i_f,i_l,nlmi,nlmj,lmj,lmi,li,lj
	pl=0.d0
	do jatom=1,natom
	  jsort=is(jatom)
	  j_f=iopb(jatom)-1
	  j_l=iopb_red(jatom)-1
	  nlmj=(lmpb_red(jsort)+1)**2
	  do iatom=1,natom
	    isort=is(iatom)
	    i_f=iopb(iatom)-1
	    i_l=iopb_red(iatom)-1
	    nlmi=(lmpb_red(isort)+1)**2
	    do lmj=1,nlmj
	      lj=lget(lmj)
	      do ie3=1,ntle_pb(lj,jsort)
	        km3=j_f+ind_prod(ie3,lmj,jsort)
	        do ie1=1,ntle_pb_red(lj,jsort)
	          km1=j_l+ind_prod_red(ie1,lmj,jsort)
	          do lmi=1,nlmi
	            li=lget(lmi)
	            do ie2=1,ntle_pb(li,isort)
	              km2=i_f+ind_prod(ie2,lmi,isort)
	              do ie=1,ntle_pb_red(li,isort)
	                km=i_l+ind_prod_red(ie,lmi,isort)
	                pl(km,km1)=pl(km,km1)+rpb_fpb(ie,ie2,li,isort)
     &				       *pf(km2,km3)*rpb_fpb(ie1,ie3,lj,jsort)
	              enddo
	            enddo
	          enddo
	        enddo
	      enddo
	    enddo
	  enddo
	enddo
      end