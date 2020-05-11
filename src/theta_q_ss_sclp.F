      subroutine theta_q_ss_sclp(k,nrow,ncol,thet_q)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nrow,ncol,k
      real*8, intent(out) :: thet_q(ncmpl*nrow,ncol,n_ixc_0)
      integer :: iatom,isort,ind,j,i,gbs_number,ix,mi,mj,iat,ii,ii0,ii1,
     &           jat,jj0,jj,jj1,jatom,jsort
      real*8 :: v(3),v1(3)
      complex*16 :: cc
      thet_q=0.d0
      do j=1,ncol
        jj=refcol_pw(j,k)
        if(jj<=n_pbmt) then
          mj=0
          do jat=1,natom
            jsort=is(jat)
            jj1=jj-iopb(jat)+1
            if(jj1>0.and.jj1<=n_pbmt0(jsort)) then
              jatom=jat
              jj0=jj1
            endif
          enddo
        else
          mj=1
          jj0=jj-n_pbmt
          jj0=indpw_gw(jj0,k)  !! G' in GBS-list
        endif
        do i=1,nrow
          ii=refrow_pw(i,k)
          if(ii<=n_pbmt) then
            mi=0
            do iat=1,natom
              isort=is(iat)
              ii1=ii-iopb(iat)+1
              if(ii1>0.and.ii1<=n_pbmt0(isort)) then
                iatom=iat
                ii0=ii1
              endif
            enddo
          else
            mi=1
            ii0=ii-n_pbmt
            ii0=indpw_gw(ii0,k)  !! G' in GBS-list
          endif
          if(mi/=mj) cycle
          if(mi==0.and.(iatom/=jatom)) cycle
          if(mi==0) then
c --------- MT part ------------------------------------------------
            do ix=1,n_ixc_0
              if(ncmpl==1) then
                thet_q(i,j,ix)=thet_me(ii0,jj0,iatom,ix)
              else
                cc=thet_me(ii0,jj0,iatom,ix)
                call zcopy(1,cc,1,thet_q(2*i-1,j,ix),1)
              endif
            enddo
          else if(mi==1) then
c --------- Interstitial -------------------------------------------
            v=gbs(:,jj0)-gbs(:,ii0)  !! G_i-G_j
            ind=gbs_number(v1)   !! G_i-G_j in GBS-list
            if(ncmpl==1) then
              do ix=1,n_ixc_0
                thet_q(i,j,ix)=theta_int(ind,ix)
              enddo
            else
              do ix=1,n_ixc_0
                call zcopy(1,theta_int(ind,ix),1,thet_q(2*i-1,j,ix),1)
              enddo
            endif
          endif
        enddo
      enddo
      end
