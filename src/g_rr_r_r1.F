      subroutine g_rr_r_r1(gk,gr)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: gk(nbasmpw,nqdiv)
      real*8, intent(out) :: gr(nqdiv,nr_full)
      integer :: k,nbas,i,i0,ind,k1,k0
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:),g0(:,:)
      allocate(tmp(nr_full))
      pi2=pi+pi
      allocate(g0(nqdiv,nr_full))
      g0=(0.d0,0.d0)
      do k=1,nqdiv
        k0=i_kref(k)
        nbas=nbask(k0)
        call zone1_number(pnt(1,k),rb0,ndiv,k1)
c ----- Fill in the first zone for first index ----------------------
        tmp=(0.d0,0.d0)
        do i=1,nbas
          i0=indgb(i,k)
          call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
          tmp(ind)=gk(i,k)
        enddo
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,tmp,1)
c ----- Multiplication with the phase factor ------------------
        do i=1,nr_full
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                             +r_pnt(3,i)*r0b0(:,3)
          fc=pi2*dot_product(pnt(:,k),v)
          cc=exp(dcmplx(0.d0,fc))
          g0(k1,i)=tmp(i)*cc
        enddo
      enddo   !! over k
      deallocate(tmp)
      call fft3(ndiv(1),ndiv(2),ndiv(3),nr_full,g0,1)
      gr=real(g0)/nqdiv
      deallocate(g0)
      end
      
      
      
      subroutine g_rr_r_r1_x(gk,gc)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: gk(nbasmpw,nrel,nqdiv)
      complex*16, intent(out) :: gc(nqdiv,nr_full,nrel,nrel)
      integer :: k,nbas,jrl,irl,ibas0,i,i0,ind,k1,k0
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:)
      allocate(tmp(nr_full))
      pi2=pi+pi
      gc=(0.d0,0.d0)
      do k=1,nqdiv
        k0=i_kref(k)
        nbas=nbask(k0)/nrel
        call zone1_number(pnt(1,k),rb0,ndiv,k1)
        do jrl=1,nrel
          do irl=1,nrel
            ibas0=(irl-1)*nbas
c ----- Fill in the first zone for first index ----------------------
            tmp=(0.d0,0.d0)
            do i=1,nbas
              i0=indgb(i,k)
              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
              tmp(ind)=gk(ibas0+i,jrl,k)
            enddo
            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,tmp,1)
c ----- Multiplication with the phase factor ------------------
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                                 +r_pnt(3,i)*r0b0(:,3)
              fc=pi2*dot_product(pnt(:,k),v)
              cc=exp(dcmplx(0.d0,fc))
              gc(k1,i,irl,jrl)=tmp(i)*cc
            enddo
          enddo   !! over irl
        enddo   !! over jrl
      enddo   !! over k
      deallocate(tmp)
      do jrl=1,nrel
        do irl=1,nrel
          call fft3(ndiv(1),ndiv(2),ndiv(3),nr_full,gc(1,1,irl,jrl),1)
        enddo
      enddo
      gc=gc/nqdiv
      end
