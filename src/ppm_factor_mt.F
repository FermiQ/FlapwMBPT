      subroutine ppm_factor_mt(ipb,fifim,z1,z2,ppm,n1,n2)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n1,n2,ipb
      real*8, intent(in) :: fifim(maxel,maxel,n_pbmtm,nsort)
      complex*16, intent(in) :: z1(nfun,n1),z2(nfun,n2)
      complex*16, intent(out) :: ppm(n1,n2)
      integer :: iatom,isort,ind,indp,ndim,ndimp,ib,i,j,km,km1
      complex*16, allocatable :: tmp(:,:)
      ppm=(0.d0,0.d0)
      allocate(tmp(maxel,n2))
      do iatom=1,natom
        isort=is(iatom)
        ind=io_lem(iatom)
        indp=iopb(iatom)-1
        ndim=lfunm(isort)
        ndimp=n_pbmt0(isort)
        do i=1,ndimp
          if(indp+i/=ipb) cycle
          tmp=(0.d0,0.d0)
          do km1=1,ndim
            j=km1+ind-1
            do km=1,ndim
              do ib=1,n2
                tmp(km,ib)=tmp(km,ib)+fifim(km,km1,i,isort)*z2(j,ib)
              enddo
            enddo
          enddo
          call zgemm('c','n',n1,n2,ndim,(1.d0,0.d0),z1(ind,1),nfun,tmp,
     &               maxel,(0.d0,0.d0),ppm,n1)
        enddo   !! over i
      enddo   !! over iatom
      deallocate(tmp)
      end
