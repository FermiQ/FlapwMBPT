      subroutine c_alpha_bands(n,k,zb,ab,ea)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: k,n
      complex*16, intent(in) :: zb(nfun,n),ab(nbasmpw,n)
      complex*16, intent(out) :: ea(n,n,3)
      integer :: nb,iatom,isort,i1,ii,iat,nt,n1,i,j,irl
      real*8 :: pi2a,v(3)
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),x(:,:,:),s(:,:),s1(:,:),
     &                           s3(:,:)
      pi2a=(pi+pi)/par
      ea=(0.d0,0.d0)
c ------- MT part --------------------------------------------------
      do isort=1,nsort
        nt=lfunm(isort)
        allocate(x(nt,nt,3))
        call c_alpha_mt(nt,x,isort)
        allocate(tmp(nt,n))
        do iat=1,nhsort(isort)
          iatom=iat_sort(iat,isort)
          i1=io_lem(iatom)
          do ii=1,3
            call zgemm('n','n',nt,n,nt,(1.d0,0.d0),x(1,1,ii),nt,
     &                 zb(i1,1),nfun,(0.d0,0.d0),tmp,nt)
            call zgemm('c','n',n,n,nt,(1.d0,0.d0),zb(i1,1),nfun,
     &                 tmp,nt,(1.d0,0.d0),ea(1,1,ii),n)
          enddo
        enddo
        deallocate(tmp,x)
      enddo
c ------- Interstitial part ----------------------------------------
      nb=nbask(k)
      if(rel_interst) then
        allocate(x(nb,nb,3))
        call c_alpha_int(nb,x,k)
        allocate(tmp(nb,n))
        do ii=1,3
          call zgemm('n','n',nb,n,nb,(1.d0,0.d0),x(1,1,ii),nb,ab(i1,1),
     &               nbasmpw,(0.d0,0.d0),tmp,nb)
          call zgemm('c','n',n,n,nb,(1.d0,0.d0),ab(1,1),nbasmpw,tmp,nb,
     &             (1.d0,0.d0),ea(1,1,ii),n)
        enddo
        deallocate(tmp,x)
      else
        n1=nb/2
        allocate(s(n1,n1))
        call s_overlap(0,s,n1,n1,indgb(1,k))
        allocate(s3(n1,n1))
        allocate(s1(n1,n))
        do j=1,3
          do i=1,n1
            v=pnt(:,k)+gbs(:,indgb(i,k))
            cc=2.d0*pi2a*v(j)/amega
            s3(:,i)=cc*s(:,i)
          enddo
          do irl=1,2
            i1=(irl-1)*n1+1
            call zgemm('n','n',n1,n,n1,(1.d0,0.d0),s3,n1,ab(i1,1),
     &                 nbasmpw,(0.d0,0.d0),s1,n1)
            call zgemm('c','n',n,n,n1,(1.d0,0.d0),ab(i1,1),nbasmpw,s1,
     &                 n1,(1.d0,0.d0),ea(1,1,j),n)
          enddo
        enddo
        deallocate(s1,s,s3)
      endif
      end
