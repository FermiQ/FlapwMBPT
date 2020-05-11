      subroutine vertex_mm_gt_r(i_tau,ll,t_rs,grs_mm,vx,fif)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll
      real*8, intent(in) :: grs_mm(nfun_red,nfun_red,0:n_tau,nqdiv_c),
     &                      fif(nindm_fif,maxel_red,maxel_red,nsort)
      complex*16, intent(in) :: t_rs(n_pbmt_red,n_pbmt_red,nqdiv_c,ll)
      complex*16, intent(inout) :: vx(nfun_red,nfun_red,nqdiv_c,ll)
      integer :: jatom,jsort,jnd,iatom,isort,ind,i,j,ii,nn,j0,n1,i0,
     &           il,nnn,j1,i1,ir,np1,n,np,il3,km,km2,iii,km1,km3
      real*8 :: tt
      real*8, allocatable :: g(:,:)
      complex*16, allocatable :: cc(:,:),aa(:,:,:),bb(:,:,:)
      do jatom=1,natom
        jsort=is(jatom)
        jnd=io_lem_red(jatom)-1
        n1=lfunm_red(jsort)
        np1=n_pbmt0_red(jsort)
	    j0=iopb_red(jatom)-1
	    do iatom=1,natom
          isort=is(iatom)
          ind=io_lem_red(iatom)-1
          n=lfunm_red(isort)
          np=n_pbmt0_red(isort)
	      i0=iopb_red(iatom)-1
	      allocate(g(n,n1))
	      allocate(aa(n,np1,n))
	      allocate(bb(n,np1,n1))
	      allocate(cc(n,n1))
          nn=n*n1
          nnn=n*np1
          do ir=1,nqdiv_c
c ----- Dynamical part : Lambda = G*T ----------------------------
            do j=1,n1
              j1=jnd+j
              do i=1,n
                i1=ind+i
                g(i,j)=grs_mm(i1,j1,i_tau,ir)
              enddo
            enddo
            do il=1,ll
              il3=3-il
              if(ll==1) il3=1
              aa=(0.d0,0.d0)
              do km=1,n
                do km2=1,n
                  do ii=1,nind_fif(km2,km,isort)
                    iii=i0+ind_fif(ii,km2,km,isort)
                    tt=fif(ii,km2,km,isort)
                    aa(km,:,km2)=aa(km,:,km2)
     &                          +tt*t_rs(j0+1:j0+np1,iii,ir,il3)
                  enddo
                enddo
              enddo
	          call dgemm('n','n',2*nnn,n1,n,1.d0,aa,2*nnn,g,n,0.d0,bb,
     &                   2*nnn)
              cc=(0.d0,0.d0)
              do km1=1,n1
                do km3=1,n1
                  do ii=1,nind_fif(km3,km1,jsort)
                    iii=ind_fif(ii,km3,km1,jsort)
                    tt=fif(ii,km3,km1,jsort)
                    cc(:,km1)=cc(:,km1)+tt*bb(:,iii,km3)
                  enddo
                enddo
              enddo
              if(ll==2) cc=cc*eqr(ir)
              do j=1,n1
                j1=jnd+j
                do i=1,n
                  i1=ind+i
                  vx(i1,j1,ir,il)=vx(i1,j1,ir,il)+cc(i,j)
                enddo
              enddo
            enddo   !! over il
          enddo   !! over ir
          deallocate(g,aa,bb,cc)
        enddo  !! over iatom
      enddo   !! over jatom
      end
