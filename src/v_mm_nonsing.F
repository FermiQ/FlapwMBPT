      subroutine v_mm_nonsing(v_full)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(inout) :: v_full(n_pbmt,n_pbmt)
      integer :: iatom,jatom,isort,jsort,ndimpb,ndimpbj,lm,l,
     &           lget,j,iei,i,indpb,jndpb,i0,lm1,l1,j0,iej
      real*8 :: a,c1,c2
      c1=16.d0*pi*pi/3.d0/amega
      c2=-32.d0*pi*pi/9.d0/amega
      do iatom=1,natom
        isort=is(iatom)
        ndimpb=(lmpb(isort)+1)**2
        indpb=iopb(iatom)-1       !! zero index in Prod basis
        do lm=1,ndimpb
          l=lget(lm)
          do iei=1,ntle_pb(l,isort)
            i0=ind_prod0(iei,l,isort)
            i=indpb+ind_prod(iei,lm,isort)
            do jatom=1,natom
              jsort=is(jatom)
              ndimpbj=(lmpb(jsort)+1)**2
              jndpb=iopb(jatom)-1       !! zero index in Prod basis
              do lm1=1,ndimpbj
                l1=lget(lm1)
                do iej=1,ntle_pb(l1,jsort)
                  j0=ind_prod0(iej,l1,jsort)
                  j=jndpb+ind_prod(iej,lm1,jsort)
                  if(l==0.and.l1==0) then
                    a=mi0(i0,isort)*mi4(iej,jsort)
     &               +mi0(j0,jsort)*mi4(iei,isort)
                    v_full(i,j)=v_full(i,j)+c1*a
                  endif
                  if(l==1.and.l1==1.and.lm==lm1) then
                    a=mi0(i0,isort)*mi0(j0,jsort)
                    v_full(i,j)=v_full(i,j)+c2*a
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      end
