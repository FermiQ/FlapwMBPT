      subroutine dos_orb(porb,weigd,ispin,ind_k,n)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
      IMPLICIT none
      integer, intent(in) :: n,ispin,ind_k
	  real*8, intent(in) :: weigd(n)
      real*8, intent(inout) :: porb(maxel_red,nsort,nspin)
      INTEGER :: k,ibnd,iatom,isort,km2,km1,ind0,jn,nn,in1,jn1,isp,k0,
     &           lf,lfm,l,ie,m,lm,k1
      real*8 :: wg
      complex*16, allocatable :: znew(:,:)
      allocate(znew(nfun,n))
      isp=min(ispin,nspin_1)
	  k0=n3_mpi_k(me_k+1)+ind_k
      do k1=1,k_star(k0)
	    k=k_list(k1,k0)
        call sym_z(znew,z_bnd(1,1,ind_k,ispin),k_group(k),pnt(1,k),nfun,
     &             n,maxb,io_lem,lmb,indbasa,limlb,0)
	    do ibnd=1,n
	      wg=weigd(ibnd)/wgt(k0)/nqdiv
          do isort=1,nsort
            iatom=iat_1(isort)
	        ind0=io_lem(iatom)-1
	        lf=0
	        lfm=0
	        do l=0,lmb(isort)
	          do nn=1,ntle_loc(l,isort)
                lf=lf+1
                ie=ind_fun0_loc(nn,l,isort)
                do m=-l,l
                  lfm=lfm+1
                  lm=l*(l+1)+m+1
			      in1=1
			      if(augm(ie,l,isort)=='LAW') in1=2
                  if(augm(ie,l,isort)=='AWL') in1=2
			      do jn1=1,in1
                    km2=ind0+indbasa(jn1,ie,lm,isort)
			        do jn=1,in1
                      km1=ind0+indbasa(jn,ie,lm,isort)
                      porb(lfm,isort,ispin)=
     &                      porb(lfm,isort,ispin)
     &		   		         +conjg(znew(km1,ibnd))*znew(km2,ibnd)
     &                       *wg*ffsmt(jn,jn1,ie,ie,l+1,isort,isp)
	                enddo
                  enddo
                enddo   !!! over m
              enddo   !!! over n
            enddo   !!! over l
          enddo        !!! over isort
        enddo               !!! over bands
      enddo
	  deallocate(znew)
      END
