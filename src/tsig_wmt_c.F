      subroutine tsig_wmt_c(g_rs_mm,k_rs,t_rs,fif,i1_tau,n,n1,np,np1,
     &                      iatom,jatom)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i1_tau,n,n1,np,np1,iatom,jatom
      complex*16, intent(in) :: fif(maxel_red,maxel_red,n_pbmtm_red,
     &                              nsort),
     &                          g_rs_mm(nfun_red,nfun_red,0:n_tau),
     &                          k_rs(n,np1)
      complex*16, intent(out) :: t_rs(np,n1)
      integer :: jsort,jnd,isort,ind,i,j,nn,ii
      complex*16 :: tt
      complex*16, allocatable :: aa(:,:,:),bb(:,:,:),dd(:,:),gg(:,:)
      isort=is(iatom)
      jsort=is(jatom)
      jnd=io_lem_red(jatom)-1
      allocate(aa(n,n1,n1))
      aa=(0.d0,0.d0)
      do j=1,n1
        do i=1,n1
          do ii=1,np1
            tt=fif(i,j,ii,jsort)
            aa(:,i,j)=aa(:,i,j)+k_rs(:,ii)*tt
          enddo
        enddo
      enddo
      ind=io_lem_red(iatom)-1
      allocate(bb(n,n1,n))
      allocate(gg(n,n1))
      gg=g_rs_mm(ind+1:ind+n,jnd+1:jnd+n1,i1_tau)
      nn=n*n1
      call zgemm('n','c',nn,n,n1,(1.d0,0.d0),aa,nn,gg,n,(0.d0,0.d0),bb,
     &           nn)
      deallocate(aa,gg)
      allocate(dd(n1,np))
      dd=(0.d0,0.d0)
      do j=1,n
        do i=1,n
          do ii=1,np
            tt=fif(i,j,ii,isort)
            dd(:,ii)=dd(:,ii)+bb(j,:,i)*tt
          enddo
        enddo
      enddo
      do j=1,n1
        do i=1,np
          t_rs(i,j)=-dd(j,i)
        enddo
      enddo
      deallocate(bb,dd)
      end
