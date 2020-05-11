      subroutine get_tmm(key,tmm)
c     key = 0 bare Coulomb interaction
c     key = 1 screened Coulomb interaction
	  use atom_mod
	  use manager_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: key
      real*8, intent(out) :: tmm(max_e_pb,max_e_pb,0:maxpb,nsort)
	  integer :: isort,l,iej,j,mtpb,iei,i,mtpbi,ir
	  real*8 :: pi8,pi8l,dqdall,c2,sqpi4
      real*8, allocatable ::  work(:),pp(:),tmp(:),f(:),g(:),bl(:,:),
     &                        cl(:,:),bc(:)
      allocate(work(0:maxnrad),pp(0:maxnrad),tmp(0:maxnrad),
     &	       f(0:maxnrad),g(0:maxnrad))
      c2=clight**2
	  pi8=8.d0*pi
      sqpi4=sqrt(pi)
	  tmm=0.d0
	  do isort=1,nsort
        if(key==1) then
          allocate(bl(0:nrad(isort),0:lmpb(isort)))
          allocate(cl(0:nrad(isort),0:lmpb(isort)))
          allocate(bc(0:lmpb(isort)))
          do ir=0,nrad(isort)
            call b_bessel(omega_hse,r(ir,isort),lmpb(isort),bc)
            bl(ir,:)=bc
            call c_bessel(omega_hse,r(ir,isort),lmpb(isort),bc)
            cl(ir,:)=bc
          enddo
          deallocate(bc)
        endif
c ------ Getting MM-values  ----------------------------------------
        do l=0,lmpb(isort)
          pi8l=pi8
	      if(key==0) pi8l=pi8/dfloat(l+l+1)
	      do iej=1,ntle_pb(l,isort)
	        j=ind_prod0(iej,l,isort)
	        mtpb=ind_pb(j,isort)
	        if(key==0) then
	          call var_int1(dr(0,isort),fun_pb(mtpb),f,g,tmp(0),
     &	                    r(0,isort),rl(0,l,isort),work,nrad(isort)+1,
     &                      h(isort))
            else
              call var_int1_b(dr(0,isort),fun_pb(mtpb),f,g,tmp(0),
     &                        r(0,isort),bl(0,l),cl(0,l),work,
     &                        nrad(isort)+1,h(isort))
            endif
	        do iei=1,ntle_pb(l,isort)
	          i=ind_prod0(iei,l,isort)
	          mtpbi=ind_pb(i,isort)
              do ir=0,nrad(isort)
                pp(ir)=tmp(ir)*fun_pb(mtpbi+ir)*dr(ir,isort)
     &		                                   *r(ir,isort)**2
              enddo
              tmm(iei,iej,l,isort)=pi8l*dqdall(h(isort),pp,nrad(isort))
	        enddo
	      enddo
	      do iej=1,ntle_pb(l,isort)
	        do iei=1,ntle_pb(l,isort)
              tmm(iei,iej,l,isort)=0.5d0*(tmm(iei,iej,l,isort)+
     &                                    tmm(iej,iei,l,isort))
              tmm(iej,iei,l,isort)=tmm(iei,iej,l,isort)
	        enddo
	      enddo
	    enddo
	    if(key==1) deallocate(bl,cl)
	  enddo   !!! over isort
      deallocate(work,pp,tmp,f,g)
      end
