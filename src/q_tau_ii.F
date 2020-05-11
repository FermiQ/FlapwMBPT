      subroutine q_tau_ii(i_tau,ll,it,il,iq,q_pw,grsii,k_rs,phf)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq,i_tau,ll,it,il
      real*8, intent(in) :: grsii(nrel,nrel*nr_red_red,nrel*nr_full_red,
     &                            0:n_tau,nqdiv_c)
      complex*16, intent(in) :: k_rs(nrel*nr_full_red,nrel*nr_full_red,
     &                               nqdiv_c,ll,2),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(inout) :: q_pw(n_pbtot_red,n_pbtot_red,nqdiv_c)
      integer :: ir,j,i,ir1,iq1,ind,kq,iq0,kq0,i1_tau,ka,i0,ig,ir0,
     &           j0,jrl,irl,js,is1,js0,is0,i1,j1,il3,ka0,kqa0,kqa,np1,
     &           np,jnd,gbs_number,id
      real*8 :: v(3),v1(3)
      complex*16 :: cc,gc
      real*8, allocatable :: sr(:,:),tr(:,:),tr1(:,:)
      complex*16, allocatable :: qt(:,:,:),t1(:,:),t2(:,:),t3(:,:),
     &                           tmp1(:,:),tmp2(:,:),s(:,:)
      i1_tau=n_tau-i_tau
      il3=ll-il+1
      allocate(qt(nqdiv_c,nr_full_red,nr_full_red))
      qt=(0.d0,0.d0)
      do i=1,nr_full_red
        i0=ir_ref_red(i)
        ig=ir_group_red(i)
        if(irel/=2) then
          do ir=1,nqdiv_c
            v=-rr_coarse(:,ir)
            call zone1_number(v,gbas,ndiv_c,ir1)
            do j=1,nr_full_red
              v=rr_coarse(:,ir)+tshift_r_red(:,i,ig)
     &                         -tshift_r_red(:,j,ig)
              call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),
     &                    u(2,ig),1)
              call zone1_number(v1,gbas,ndiv_c,ir0)
              j0=ip_r_red(j,ig)
c ------------------------ A*B part ------------------------------
              qt(ir,i,j)=qt(ir,i,j)+grsii(1,i0,j0,i1_tau,ir0)
     &                             *k_rs(i,j,ir,il,it)
c ------------------------ C*D part ------------------------------
              cc=grsii(1,i0,j0,i_tau,ir0)*k_rs(j,i,ir1,il3,3-it)
              if(ll==2) cc=cc*eqr(ir)
              qt(ir,i,j)=qt(ir,i,j)+cc
            enddo  !! over j
          enddo   !! over ir
        else
          do ir=1,nqdiv_c
            do j=1,nr_full_red 
              v=rr_coarse(:,ir)+tshift_r_red(:,i,ig)
     &                         -tshift_r_red(:,j,ig)
              call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
              call zone1_number(v1,gbas,ndiv_c,ir0)
              j0=ip_r_red(j,ig)
              do jrl=1,2
                j1=(jrl-1)*nr_full_red+j
                do irl=1,2
                  i1=(irl-1)*nr_full_red+i
                  do js=1,2
                    js0=(js-1)*nr_full_red+j0
                    do is1=1,2
                      is0=(is1-1)*nr_red_red+i0
                      cc=conjg(uj((is1-1)*2+irl,ig))
     &                        *uj((js-1)*2+jrl,ig)
                      if(ll==2) cc=cc*eqr(ir)
c ------------------------ A*B part ------------------------------
                      gc=dcmplx(grsii(1,is0,js0,i1_tau,ir0),
     &                         -grsii(2,is0,js0,i1_tau,ir0))
                      qt(ir,i,j)=qt(ir,i,j)+cc*gc
     &                                     *k_rs(i1,j1,ir,il,it)
c ------------------------ C*D part ------------------------------
                      gc=dcmplx(grsii(1,is0,js0,i_tau,ir0),
     &                          grsii(2,is0,js0,i_tau,ir0))
                      qt(ir,i,j)=qt(ir,i,j)+cc*gc
     &                                     *k_rs(i1,j1,ir1,il3,3-it)
                    enddo
                  enddo
                enddo
              enddo
            enddo   !! over j
          enddo   !! over ir
        endif
      enddo   !! over i
      call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nr_full_red**2,qt,-1)
      allocate(t1(nr_full_red,nr_full_red))
      allocate(t2(nr_full_red,nplw_gw_red))
      do iq1=1,nqdiv_c
        iq0=i_kref_c(iq1)
        ka0=k_a_from_c(iq0)
        ka=k_a_from_c(iq1)
        np=nplwgw_red(ka0)
        v=pnt_c(:,iq1)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)     !! for Q'-Q
        kq0=i_kref_c(kq)
        kqa0=k_a_from_c(kq0)
        kqa=k_a_from_c(kq)
        np1=nplwgw_red(kqa0)
        call zone1_number(pnt_c(1,iq1),rb0_c,ndiv_c,ind)
        do ir=1,nr_full_red
          cc=phf(ir,kq)
          t1(ir,:)=cc*qt(ind,:,ir)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nr_full_red,t1,
     &            1)
        do i=1,np1
          i0=indpw_gw_red(i,kqa)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          t2(:,i)=t1(ind,:)
        enddo
        t2=t2/nr_full_red*sqrt(amega)
        do ir=1,nr_full_red
          cc=conjg(phf(ir,iq1))
          t2(ir,1:np1)=cc*t2(ir,1:np1)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),np1,t2,-1)
        allocate(t3(np,np1))
        do i=1,np
          i0=indpw_gw_red(i,ka)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          t3(i,1:np1)=t2(ind,1:np1)
        enddo
        t3=t3/nr_full_red*sqrt(amega)
        allocate(tmp2(np,np1))
        if(complex_ro) then
          allocate(s(np,np))
          do j=1,np
            jnd=indpw_gw_red(j,ka)
            do i=1,np
              ind=indpw_gw_red(i,ka)
              v=gbs(:,jnd)-gbs(:,ind)
              id=gbs_number(v)
              s(i,j)=dcmplx(sovr(id),sovi(id))
            enddo
          enddo
          allocate(tmp1(np,np1))
          call zgemm('n','n',np,np1,np,(1.d0,0.d0),s,np,t3,np,
     &               (0.d0,0.d0),tmp1,np)
          deallocate(s)
          allocate(s(np1,np1))
          do j=1,np1
            jnd=indpw_gw_red(j,kqa)
            do i=1,np1
              ind=indpw_gw_red(i,kqa)
              v=gbs(:,jnd)-gbs(:,ind)
              id=gbs_number(v)
              s(i,j)=dcmplx(sovr(id),sovi(id))
            enddo
          enddo
          call zgemm('n','n',np,np1,np1,(1.d0,0.d0),tmp1,np,s,np1,
     &               (0.d0,0.d0),tmp2,np)
          deallocate(s,tmp1)
        else
          allocate(sr(np1,np1))
          do j=1,np1
            jnd=indpw_gw_red(j,kqa)
            do i=1,np1
              ind=indpw_gw_red(i,kqa)
              v=gbs(:,jnd)-gbs(:,ind)
              id=gbs_number(v)
              sr(i,j)=sovr(id)
            enddo
          enddo
          allocate(tmp1(np,np1))
          call dgemm('n','n',2*np,np1,np1,1.d0,t3,2*np,sr,np1,0.d0,
     &               tmp1,2*np)
          deallocate(sr)
          allocate(sr(np,np))
          do j=1,np
            jnd=indpw_gw_red(j,ka)
            do i=1,np
              ind=indpw_gw_red(i,ka)
              v=gbs(:,jnd)-gbs(:,ind)
              id=gbs_number(v)
              sr(i,j)=sovr(id)
            enddo
          enddo
          allocate(tr(np,np1))
          allocate(tr1(np,np1))
          tr=real(tmp1)
          call dgemm('n','n',np,np1,np,1.d0,sr,np,tr,np,0.d0,tr1,np)
          tmp2=tr1
          tr=imag(tmp1)
          call dgemm('n','n',np,np1,np,1.d0,sr,np,tr,np,0.d0,tr1,np)
          tmp2=tmp2+(0.d0,1.d0)*tr1
          deallocate(sr,tr,tr1,tmp1)
        endif
        deallocate(t3)
        do j=1,np1
          do i=1,np
            q_pw(n_pbmt_red+i,n_pbmt_red+j,iq1)=-tmp2(i,j)
          enddo
        enddo
        deallocate(tmp2)
      enddo   !! over iq1
      deallocate(qt,t1,t2)
      end
