	subroutine sym_left_red(ig,f0,f)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
	use vertex_mod
	implicit none
	integer, intent(in) :: ig
	complex*16, intent(in) :: f0(nfun_red)
	complex*16, intent(out) :: f(nfun_red)
	integer :: iatom,jatom,isort,ind0,jnd0,l,km,lm0,lf,li,ii,in,
     &	         ie,lf0,jj,mj,lm
	complex*16 :: y0(2*maxb_red+nrel),t0(2*maxb_red+nrel)
	do iatom=1,natom
	  jatom=ip(iatom,ig)
	  isort=is(iatom)
	  ind0=io_lem_red(iatom)-1   !! zero index for atom iatom_c
	  jnd0=io_lem_red(jatom)-1
	  do lf=1,lfun_red(isort)
	    lf0=le_red(lf,isort)
	    li=li_val(lf0,isort)
	    ie=ie_val(lf0,isort)
	    in=in_val(lf0,isort)
	    if(irel<=1) then
	      l=li-1
	      jj=l
	    else
		    l=li/2
		    ii=li-2*l
	      if(ii.eq.0) ii=-1
            jj=l+l+ii
          endif
          lm0=0
          do mj=-jj,jj,nrel
            lm0=lm0+1
            if(irel.lt.2) then
			  lm=l*(l+1)+mj+1
	      else
		      call getlimj(lm,l,ii,mj,li,1)
	      endif
            km=indbasa_red(in,ie,lm,isort)
	      t0(lm0)=f0(jnd0+km)
	    enddo
	    if(irel<=1) call rotate1c(t0,y0,l,u(1,ig),1)
	    if(irel==2) call rotate2(t0,y0,l,ii,uj(1,ig),1,1)
          lm0=0
          do mj=-jj,jj,nrel
            lm0=lm0+1
            if(irel.lt.2) then
			  lm=l*(l+1)+mj+1
	      else
		      call getlimj(lm,l,ii,mj,li,1)
	      endif
            km=indbasa_red(in,ie,lm,isort)
	      f(ind0+km)=y0(lm0)
	    enddo
	  enddo   !! over l
	enddo   !! over iatom
	end