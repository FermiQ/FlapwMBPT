      subroutine q_tau_mm_c(i_tau,ll,it,il,q_pw,grsmm,k_rs,fif,n0p)
	use atom_mod
	use manager_mod
      use models_mod
	use parallel_mod
	use solid_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll,it,il,n0p
      complex*16, intent(in) :: k_rs(nfun_red,nfun_red,nqdiv_c,ll,2),
     &                          grsmm(nfun_red,nfun_red,0:n_tau,
     &                                nqdiv_c),
     &                          fif(nindm_fif,maxel_red,maxel_red,nsort)
      complex*16, intent(inout) :: q_pw(n0p,n0p,nqdiv_c)
      integer :: ir,j,i,ir1,iq1,jatom,jsort,iatom,isort,j0,jnd,n1,i0,
     &           ind,n,np1,np,iqq,i1_tau,il3,km1,km3,ii,iii,km,km2
      real*8 :: v(3)
      complex*16 :: tt
      complex*16, allocatable :: qt(:,:,:),aa(:,:,:),tmp(:,:),bb(:,:,:)
      i1_tau=n_tau-i_tau
      il3=ll-il+1
	do jatom=1,natom
        jsort=is(jatom)
        j0=iopb_red(jatom)-1
        jnd=io_lem_red(jatom)-1
        n1=lfunm_red(jsort)
        np1=n_pbmt0_red(jsort)
        do iatom=1,natom
          isort=is(iatom)
          i0=iopb_red(iatom)-1
          ind=io_lem_red(iatom)-1
          n=lfunm_red(isort)
          np=n_pbmt0_red(isort)
          allocate(aa(np1,n,n1))
          allocate(bb(n,n1,np))
          allocate(qt(nqdiv_c,np,np1))
          allocate(tmp(np1,np))
          qt=(0.d0,0.d0)
	    do ir=1,nqdiv_c
            v=-rr_coarse(:,ir)
            call zone1_number(v,gbas,ndiv_c,ir1)
c --------------------- A*B part ------------------------------
            aa=(0.d0,0.d0)
            do km1=1,n1
              do km3=1,n1
                do ii=1,nind_fif(km3,km1,jsort)
                  iii=ind_fif(ii,km3,km1,jsort)
                  tt=fif(ii,km3,km1,jsort)
                  aa(iii,:,km1)=aa(iii,:,km1)+tt*k_rs(ind+1:ind+n,
     &                                                jnd+km3,ir,il,it)
                enddo
              enddo
            enddo
            bb=(0.d0,0.d0)
            do km=1,n
              do km2=1,n
                do ii=1,nind_fif(km2,km,isort)
                  iii=ind_fif(ii,km2,km,isort)
                  tt=fif(ii,km2,km,isort)
                  bb(km2,:,iii)=bb(km2,:,iii)
     &                         +tt*grsmm(ind+km,jnd+1:jnd+n1,i1_tau,ir)
                enddo
              enddo
            enddo
            bb=conjg(bb)
	      call zgemm('n','n',np1,np,n*n1,(1.d0,0.d0),aa,np1,bb,n*n1,
     &                 (0.d0,0.d0),tmp,np1)
            do j=1,np1
              do i=1,np
                qt(ir,i,j)=qt(ir,i,j)+tmp(j,i)
              enddo
            enddo
c ------------------------ C*D part ------------------------------
            aa=(0.d0,0.d0)
            do km1=1,n1
              do km3=1,n1
                do ii=1,nind_fif(km3,km1,jsort)
                  iii=ind_fif(ii,km3,km1,jsort)
                  tt=fif(ii,km3,km1,jsort)
                  aa(iii,:,km1)=aa(iii,:,km1)
     &                    +tt*k_rs(jnd+km3,ind+1:ind+n,ir1,il3,3-it)
                enddo
              enddo
            enddo
            if(ll==2) aa=aa*eqr(ir)
            bb=(0.d0,0.d0)
            do km=1,n
              do km2=1,n
                do ii=1,nind_fif(km2,km,isort)
                  iii=ind_fif(ii,km2,km,isort)
                  tt=conjg(fif(ii,km2,km,isort))
                  bb(km2,:,iii)=bb(km2,:,iii)
     &                         +tt*grsmm(ind+km,jnd+1:jnd+n1,i_tau,ir)
                enddo
              enddo
            enddo
	      call zgemm('n','n',np1,np,n*n1,(1.d0,0.d0),aa,np1,bb,n*n1,
     &                 (0.d0,0.d0),tmp,np1)
            do j=1,np1
              do i=1,np
                qt(ir,i,j)=qt(ir,i,j)+tmp(j,i)
              enddo
            enddo
          enddo   !! over ir
          if(n0p==n_pbtot_red) then
            call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),np*np1,qt,-1)
            do iq1=1,nqdiv_c
              call zone1_number(pnt_c(1,iq1),rb0_c,ndiv_c,iqq)
              do j=1,np1
                do i=1,np
                  q_pw(i0+i,j0+j,iq1)=-qt(iqq,i,j)
                enddo
              enddo
            enddo   !! over iq1
          else if(n0p==n_pbmt_red) then
            do ir=1,nqdiv_c
              do j=1,np1
                do i=1,np
                  q_pw(i0+i,j0+j,ir)=-qt(ir,i,j)
                enddo
              enddo
            enddo   !! over iq1
          endif
          deallocate(aa,qt,tmp,bb)
        enddo   !! over iatom_c
      enddo  !! over jatom_c
      end