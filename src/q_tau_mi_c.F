      subroutine q_tau_mi_c(i_tau,ll,it,il,iq,q_pw,grsmi,grsim,k_rs,fif,
     &                      phf)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,iq,ll,it,il
      complex*16, intent(in) :: grsmi(nfun_red,2*nr_full_red,0:n_tau,
     &                                nqdiv_c),
     &                          grsim(2*nr_full_red,nfun_red,0:n_tau,
     &                                nqdiv_c),
     &                          fif(nindm_fif,maxel_red,maxel_red,
     &                              nsort),
     &                          k_rs(nfun_red,2*nr_full_red,nqdiv_c,ll,
     &                               2),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(inout) :: q_pw(n_pbtot_red,n_pbtot_red,nqdiv_c)
      integer :: ir,j,i,ir1,iq1,iatom,isort,i0,ind,n,np,is1,iq0,iqq,
     &           i1_tau,ka,irl,jj,kq,il3,km1,km3,ii,iii,km,km2,ka0,npg,
     &           jnd,id,gbs_number,kq0,ip0
      real*8 :: v(3)
      complex*16 :: cc,tt
      real*8, allocatable :: sr(:,:),tr(:,:),tr1(:,:)
      complex*16, allocatable :: qt(:,:,:),aa(:),t1(:,:),t2(:,:),
     &                           qt1(:,:,:),s(:,:),tmp1(:,:),bb(:,:,:),
     &                           bb1(:,:,:)
      i1_tau=n_tau-i_tau
      il3=ll-il+1
      do iatom=1,natom
        isort=is(iatom)
        ip0=iopb_red(iatom)-1
        ind=io_lem_red(iatom)-1
        n=lfunm_red(isort)
        np=n_pbmt0_red(isort)
        allocate(aa(n))
        allocate(bb1(2*nr_full_red,n,np))
        allocate(bb(n,2*nr_full_red,np))
        allocate(qt(nqdiv_c,np,nr_full_red))
        allocate(qt1(nqdiv_c,nr_full_red,np))
        qt1=(0.d0,0.d0)
        qt=(0.d0,0.d0)
        do ir=1,nqdiv_c
          v=-rr_coarse(:,ir)
          call zone1_number(v,gbas,ndiv_c,ir1)
c -------------- Contribution to Q_im ------------------------------
          bb1=(0.d0,0.d0)
          do km1=1,n
            do km3=1,n
              do ii=1,nind_fif(km3,km1,isort)
                iii=ind_fif(ii,km3,km1,isort)
                tt=fif(ii,km3,km1,isort)
                bb1(:,km3,iii)=bb1(:,km3,iii)
     &                        +tt*conjg(grsim(:,ind+km1,i_tau,ir))
              enddo
            enddo
          enddo  
          do j=1,nr_full_red
            do irl=1,2
              jj=(irl-1)*nr_full_red+j
              aa=k_rs(ind+1:ind+n,jj,ir1,il3,3-it)
              if(ll==2) aa=aa*eqr(ir)
              do is1=1,np
                cc=dot_product(bb1(jj,:,is1),aa)
                qt1(ir,j,is1)=qt1(ir,j,is1)+cc
              enddo
            enddo
          enddo         
c ---------- Contribution to Q_mi ------------------------------
          bb=(0.d0,0.d0)
          do km=1,n
            do km2=1,n
              do ii=1,nind_fif(km2,km,isort)
                iii=ind_fif(ii,km2,km,isort)
                tt=fif(ii,km2,km,isort)
                bb(km2,:,iii)=bb(km2,:,iii)
     &                       +tt*conjg(grsmi(ind+km,:,i1_tau,ir))
              enddo
            enddo
          enddo   
          do j=1,nr_full_red
            do irl=1,2
              jj=(irl-1)*nr_full_red+j
              aa=k_rs(ind+1:ind+n,jj,ir,il,it)
              do is1=1,np
                cc=dot_product(bb(:,jj,is1),aa)
                qt(ir,is1,j)=qt(ir,is1,j)+cc
              enddo
            enddo
          enddo
        enddo   !! over ir
        call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),np*nr_full_red,qt,-1)
        call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),np*nr_full_red,qt1,-1)
        deallocate(aa,bb,bb1)
        allocate(t1(nr_full_red,np))
        do iq1=1,nqdiv_c
          iq0=i_kref_c(iq1)
          ka=k_a_from_c(iq1)
          ka0=k_a_from_c(iq0)
          call zone1_number(pnt_c(1,iq1),rb0_c,ndiv_c,iqq)
          do ir=1,nr_full_red
            cc=conjg(phf(ir,iq1))
            t1(ir,:)=cc*qt1(iqq,ir,:)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),np,t1,-1)
          npg=nplwgw_red(ka0)
          allocate(t2(npg,np))
          t2=(0.d0,0.d0)
          do i=1,npg
            i0=indpw_gw_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t2(i,:)=t1(ind,:)
          enddo
          t2=t2/nr_full_red*sqrt(amega)
          allocate(tmp1(npg,np))
          if(complex_ro) then
            allocate(s(npg,npg))
            do j=1,npg
              jnd=indpw_gw_red(j,ka)
              do i=1,npg
                ind=indpw_gw_red(i,ka)
                v=gbs(:,jnd)-gbs(:,ind)
                id=gbs_number(v)
                s(i,j)=dcmplx(sovr(id),sovi(id))
              enddo
            enddo
            call zgemm('n','n',npg,np,npg,(1.d0,0.d0),s,npg,t2,npg,
     &                 (0.d0,0.d0),tmp1,npg)
            deallocate(s)
          else
            allocate(sr(npg,npg))
            do j=1,npg
              jnd=indpw_gw_red(j,ka)
              do i=1,npg
                ind=indpw_gw_red(i,ka)
                v=gbs(:,jnd)-gbs(:,ind)
                id=gbs_number(v)
                sr(i,j)=sovr(id)
              enddo
            enddo
            allocate(tr(npg,np))
            allocate(tr1(npg,np))
            tr=real(t2)
            call dgemm('n','n',npg,np,npg,1.d0,sr,npg,tr,npg,0.d0,tr1,
     &                 npg)
            tmp1=tr1
            tr=imag(t2)
            call dgemm('n','n',npg,np,npg,1.d0,sr,npg,tr,npg,0.d0,tr1,
     &                 npg)
            tmp1=tmp1+(0.d0,1.d0)*tr1
            deallocate(sr,tr,tr1)
          endif
          deallocate(t2)
          do j=1,np
            do i=1,npg
              q_pw(n_pbmt_red+i,ip0+j,iq1)=q_pw(n_pbmt_red+i,ip0+j,iq1)
     &                                    -tmp1(i,j)
            enddo
          enddo
          deallocate(tmp1)
        enddo   !! over iq1
        deallocate(t1)
        allocate(t1(nr_full_red,np))
        do iq1=1,nqdiv_c
          v=pnt_c(:,iq1)-pnt_c(:,iq)
          call zone1_number(v,rb0_c,ndiv_c,kq)
          kq=index_k1_c(kq)     !! for Q'-Q
          kq0=i_kref_c(kq)
          ka0=k_a_from_c(kq0)
          ka=k_a_from_c(kq)
          call zone1_number(pnt_c(1,iq1),rb0_c,ndiv_c,iqq)
          do ir=1,nr_full_red
            cc=phf(ir,kq)
            t1(ir,:)=cc*qt(iqq,:,ir)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),np,t1,1)
          npg=nplwgw_red(ka0)
          allocate(t2(np,npg))
          t2=(0.d0,0.d0)
          do i=1,npg
            i0=indpw_gw_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t2(:,i)=t1(ind,:)
          enddo
          t2=t2/nr_full_red*sqrt(amega)
          allocate(tmp1(np,npg))
          if(complex_ro) then
            allocate(s(npg,npg))
            do j=1,npg
              jnd=indpw_gw_red(j,ka)
              do i=1,npg
                ind=indpw_gw_red(i,ka)
                v=gbs(:,jnd)-gbs(:,ind)
                id=gbs_number(v)
                s(i,j)=dcmplx(sovr(id),sovi(id))
              enddo
            enddo
            call zgemm('n','n',np,npg,npg,(1.d0,0.d0),t2,np,s,npg,
     &                 (0.d0,0.d0),tmp1,np)
            deallocate(s)
          else
            allocate(sr(npg,npg))
            do j=1,npg
              jnd=indpw_gw_red(j,ka)
              do i=1,npg
                ind=indpw_gw_red(i,ka)
                v=gbs(:,jnd)-gbs(:,ind)
                id=gbs_number(v)
                sr(i,j)=sovr(id)
              enddo
            enddo
            allocate(tr(np,npg))
            allocate(tr1(np,npg))
            tr=real(t2)
            call dgemm('n','n',np,npg,npg,1.d0,tr,np,sr,npg,0.d0,tr1,
     &                 np)
            tmp1=tr1
            tr=imag(t2)
            call dgemm('n','n',np,npg,npg,1.d0,tr,np,sr,npg,0.d0,tr1,
     &                 npg)
            tmp1=tmp1+(0.d0,1.d0)*tr1
            deallocate(sr,tr,tr1)
          endif
          deallocate(t2)
          do j=1,npg
            do i=1,np
              q_pw(ip0+i,n_pbmt_red+j,iq1)=q_pw(ip0+i,n_pbmt_red+j,iq1)
     &                                    -tmp1(i,j)
            enddo
          enddo
          deallocate(tmp1)
        enddo   !! over iq1
        deallocate(qt,t1,qt1)
      enddo   !! over iatom_c 
      end
