      subroutine sym_w_mt_sort(key,ig,tr,w_0,n,nat,npbmt,isort)
c     See identity VII-32-3
c     key = 1 - left transformation
c     key = 2 - right transformation
      use atom_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ig,key,n,nat,npbmt,isort
      complex*16, intent(in) :: tr(nat)
      complex*16, intent(inout) :: w_0(npbmt,n)
      integer :: i,iatom,i0
      real*8, allocatable :: sr(:,:),sj(:,:),s0(:,:)
      allocate(s0(npbmt,n),sj(npbmt,n),sr(npbmt,n))
      s0=real(w_0)
      call sym_pro_w_sort(s0,sr,ip(1,ig),u(1,ig),n,npbmt,isort,nat)
      s0=imag(w_0)
      call sym_pro_w_sort(s0,sj,ip(1,ig),u(1,ig),n,npbmt,isort,nat)
      w_0=dcmplx(sr,sj)
      deallocate(s0,sr,sj)
      do iatom=1,nat
        i0=n_pbmt0(isort)*(iatom-1)+1
        i=i0+n_pbmt0(isort)-1
        if(key==1) w_0(i0:i,:)=tr(iatom)*w_0(i0:i,:)
        if(key==2) w_0(i0:i,:)=conjg(tr(iatom))*w_0(i0:i,:)
      enddo
      end
