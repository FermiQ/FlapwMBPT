      subroutine v_coul_new_mt(v_full,vmm,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: khbr
      complex*16, intent(in) :: vmm(nlbs,nlbs)
      complex*16, intent(out) :: v_full(n_pbmt,n_pbmt)
      integer :: iatom,jatom,isort,jsort,ndimpb,ndimpbj,ind,jnd,lm,l,
     &           lget,j,iei,i,indpb,jndpb,indi,indj,i0,lm1,l1,j0,iej
c --------------------------------------------------------------------
      v_full=(0.d0,0.d0)
c ------ Filling the local part of V_TMP ------------------------------
      do iatom=1,natom
        isort=is(iatom)
        ndimpb=(lmpb(isort)+1)**2
        ind=iopb(iatom)-1
        do lm=1,ndimpb
          l=lget(lm)
          do iej=1,ntle_pb(l,isort)
            j=ind+ind_prod(iej,lm,isort)
            do iei=1,ntle_pb(l,isort)
              i=ind+ind_prod(iei,lm,isort)
              if(khbr==0) v_full(i,j)=tmm_pb(iei,iej,l,isort)
              if(khbr==1) v_full(i,j)=tmm_pb_b(iei,iej,l,isort)
            enddo
          enddo
        enddo
      enddo
c ------ Adding the nonlocal MT-MT part to V_TMP -------------------
      do iatom=1,natom
        isort=is(iatom)
        ndimpb=(lmpb(isort)+1)**2
        ind=iobs(iatom)-1       !! zero index in V basis
        indpb=iopb(iatom)-1       !! zero index in Prod basis
        do lm=1,ndimpb
          l=lget(lm)
          indi=ind+lm
          do iei=1,ntle_pb(l,isort)
            i0=ind_prod0(iei,l,isort)
            i=indpb+ind_prod(iei,lm,isort)
            do jatom=1,natom
              jsort=is(jatom)
              ndimpbj=(lmpb(jsort)+1)**2
              jnd=iobs(jatom)-1       !! zero index in V basis
              jndpb=iopb(jatom)-1       !! zero index in Prod basis
              do lm1=1,ndimpbj
                l1=lget(lm1)
                indj=jnd+lm1
                do iej=1,ntle_pb(l1,jsort)
                  j0=ind_prod0(iej,l1,jsort)
                  j=jndpb+ind_prod(iej,lm1,jsort)
                  if(khbr==0) v_full(i,j)=v_full(i,j)+mi0(i0,isort)
     &                                *vmm(indi,indj)*mi0(j0,jsort)
                  if(khbr==1) v_full(i,j)=v_full(i,j)+mi0_b(i0,isort)
     &                                *vmm(indi,indj)*mi0_b(j0,jsort)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      end
