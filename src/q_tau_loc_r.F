      subroutine q_tau_loc_r(i_tau,ll,it,il,q_pw,gtau,k_rs,fif,n,np)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll,it,il,n,np
      real*8, intent(in) :: gtau(n,n,0:n_tau),fif(n,n,np)
      complex*16, intent(in) :: k_rs(n,n,ll,2)
      complex*16, intent(inout) :: q_pw(np,np)
      integer :: i1_tau,il3,ii
      real*8, allocatable :: aa(:,:,:),b(:,:,:),c(:,:)
      complex*16, allocatable :: dd(:,:),bb(:,:,:),tmp(:,:)
      i1_tau=n_tau-i_tau
      il3=ll-il+1
      allocate(aa(n,n,np))
      allocate(bb(n,n,np))
      allocate(tmp(np,np))
      allocate(b(n,n,np))
      allocate(c(np,np))
      allocate(dd(n,n))
c --------------------- A*B part ------------------------------
      do ii=1,np
        call dgemm('n','n',n,n,n,1.d0,fif(1,1,ii),n,gtau(1,1,i1_tau),n,
     &             0.d0,aa(1,1,ii),n)
      enddo
      call dgemm('n','n',2*n,n*np,n,1.d0,k_rs(1,1,il,it),2*n,fif,n,0.d0,
     &           bb,2*n)
      b=real(bb)
      call dgemm('t','n',np,np,n*n,1.d0,aa,n*n,b,n*n,0.d0,c,np)
      tmp=c
      b=imag(bb)
      call dgemm('t','n',np,np,n*n,1.d0,aa,n*n,b,n*n,0.d0,c,np)
      tmp=tmp+(0.d0,1.d0)*c
      q_pw=q_pw-tmp
c ------------------------ C*D part ------------------------------
      do ii=1,np
        call dgemm('n','n',n,n,n,1.d0,fif(1,1,ii),n,gtau(1,1,i_tau),n,
     &             0.d0,aa(1,1,ii),n)
      enddo
      dd=transpose(k_rs(:,:,il3,3-it))
      call dgemm('n','n',2*n,n*np,n,1.d0,dd,2*n,fif,n,0.d0,bb,2*n)
      b=real(bb)
      call dgemm('t','n',np,np,n*n,1.d0,aa,n*n,b,n*n,0.d0,c,np)
      tmp=c
      b=imag(bb)
      call dgemm('t','n',np,np,n*n,1.d0,aa,n*n,b,n*n,0.d0,c,np)
      tmp=tmp+(0.d0,1.d0)*c
      q_pw=q_pw-tmp
      deallocate(aa,bb,tmp,b,c,dd)
      end
