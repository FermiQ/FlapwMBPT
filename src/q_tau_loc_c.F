      subroutine q_tau_loc_c(i_tau,ll,it,il,q_pw,gtau,k_rs,fif,n,np)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll,it,il,n,np
      complex*16, intent(in) :: gtau(n,n,0:n_tau),fif(n,n,np),
     &                          k_rs(n,n,ll,2)
      complex*16, intent(inout) :: q_pw(np,np)
      integer :: i1_tau,il3,ii
      complex*16, allocatable :: aa(:,:,:),dd(:,:),bb(:,:,:),tmp(:,:)
      i1_tau=n_tau-i_tau
      il3=ll-il+1
      allocate(aa(n,n,np))
      allocate(bb(n,n,np))
      allocate(tmp(np,np))
      allocate(dd(n,n))
c --------------------- A*B part ------------------------------
      do ii=1,np
        call zgemm('n','n',n,n,n,(1.d0,0.d0),fif(1,1,ii),n,
     &             gtau(1,1,i1_tau),n,(0.d0,0.d0),aa(1,1,ii),n)
      enddo
      aa=conjg(aa)
      call zgemm('n','n',n,n*np,n,(1.d0,0.d0),k_rs(1,1,il,it),n,fif,n,
     &           (0.d0,0.d0),bb,n)
      call zgemm('c','n',np,np,n*n,(1.d0,0.d0),aa,n*n,bb,n*n,
     &           (0.d0,0.d0),tmp,np)
      q_pw=q_pw-tmp
c ------------------------ C*D part ------------------------------
      bb=conjg(fif)
      do ii=1,np
        call zgemm('n','n',n,n,n,(1.d0,0.d0),bb(1,1,ii),n,
     &             gtau(1,1,i_tau),n,(0.d0,0.d0),aa(1,1,ii),n)
      enddo
      dd=transpose(k_rs(:,:,il3,3-it))
      call zgemm('n','n',n,n*np,n,(1.d0,0.d0),dd,n,fif,n,(0.d0,0.d0),bb,
     &           n)
      call zgemm('c','n',np,np,n*n,(1.d0,0.d0),aa,n*n,bb,n*n,
     &           (0.d0,0.d0),tmp,np)
      q_pw=q_pw-tmp
      q_pw=q_pw-tmp
      deallocate(aa,bb,tmp,dd)
      end
