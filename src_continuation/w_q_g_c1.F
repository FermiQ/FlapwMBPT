      subroutine w_q_g_c1(ind_tau,w_c)
C Calculates W(q;G;r') where q=1,npnt+par; G=1,nplw_gw; r'=1,nr_full -
      use atom_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau
      complex*16, intent(out) :: w_c(nplw_gw,npnt,max_r_star,
     &                               ndim4_nr_red(me_k_b+1))
      integer :: k,ind_k,j,j0,i,ind,ind_ista,ista,npw,n0,n1,n2,i0,ij,
     &           ind_i,nn,ii,jj,ij0
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      integer, allocatable :: nd_b(:),n0_b(:)
      complex*16, allocatable :: tmp(:,:),w(:),wc(:,:,:)
      pi2=pi+pi
      allocate(nd_b(nproc_b))
      allocate(n0_b(nproc_b))
      call size_shift_par(nplw_gw,nproc_b,nd_b,n0_b)
      n0=n0_b(me_b+1)
      n1=n0+nd_b(me_b+1)
      allocate(wc(nd_b(me_b+1),ndim3_k(me_k+1),nr_full))
      wc=(0.d0,0.d0)
      allocate(tmp(nr_full,nd_b(me_b+1)))
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        npw=nplwgw(k)
        if(npw==0) cycle
        nn=npw*(npw+1)/2
        allocate(w(nn))
        w=(0.d0,0.d0)
        if(lsclp) then
          do j=1,ncol_pw(k)
            jj=refcol_pw(j,k)-n_pbmt
            if(jj<1) cycle
            ij0=jj*(jj-1)/2
            do i=1,nrow_pw(k)
              ii=refrow_pw(i,k)-n_pbmt
              if(ii<1.or.ii>jj) cycle
              ij=ij0+ii
              call zcopy(1,pw_sclp(2*i-1,j,ind_tau,ind_k),1,w(ij),1)
            enddo
          enddo
        else
          do ind_i=1,nd_b_pbint(me_b+1)
            i=n0_b_pbint(me_b+1)+ind_i
            if(i<=nn) call zcopy(1,
     &                           p_w_q_tau_ii(2*ind_i-1,ind_tau,ind_k),
     &                           1,w(i),1)
          enddo
        endif
        if(nproc_b/=1) call dgop(w,2*nn,'  +',comm_b)
        if(n0>=npw) then
          deallocate(w)
          cycle
        endif
        n2=min(nd_b(me_b+1),npw-n0)
        w=w/2
        tmp=(0.d0,0.d0)
        do j=1,npw
          j0=indpw_gw(j,k)
          call zone1_number(gbs(:,j0),rbas,nrdiv,ind)
          do i=1,n2
            i0=n0+i
            if(i0<=j) then
              ij=j*(j-1)/2+i0
              tmp(ind,i)=w(ij)
            else
              ij=i0*(i0-1)/2+j
              tmp(ind,i)=conjg(w(ij))
            endif
          enddo
        enddo
        deallocate(w)
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n2,tmp,-1)
c ----- Multiplication with the phase factor ------------------
        do i=1,nr_full
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                             +r_pnt(3,i)*r0b0(:,3)
          fc=-pi2*dot_product(pnt(:,k),v)
          cc=exp(dcmplx(0.d0,fc))/sqrt(amega)
          wc(1:n2,ind_k,i)=tmp(i,1:n2)*cc
        enddo
      enddo   !! over ind_k
      deallocate(tmp)
      allocate(tmp(npnt,nr_full))
      do i=1,nplw_gw
        tmp=(0.d0,0.d0)
        if(i>n0.and.i<=n1) then
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            tmp(k,:)=wc(i-n0,ind_k,:)
          enddo
        endif
        if(nproc_k*nproc_b/=1) call dgop(tmp,2*nr_full*npnt,'  +',
     &                                   comm_k_b)
        do ind_ista=1,ndim4_nr_red(me_k_b+1)
          ista=n4_mpi_nr_red(me_k_b+1)+ind_ista
          do j=1,num_r_star(ista)
            ind=ind_r_star(j,ista)
            w_c(i,:,j,ind_ista)=tmp(:,ind)
          enddo
        enddo
      enddo
      deallocate(wc,tmp,nd_b,n0_b)
      end
