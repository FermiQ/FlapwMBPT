      subroutine mt_integral_neighbors(mtp,pb,iq,chi_neig)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq
      real*8, intent(in) :: mtp(n_pbmt)
      complex*16, intent(in) :: pb(n_pbmt,n_pbmt)
      complex*16, intent(inout) :: chi_neig(max_sites,nsort)
      integer :: iatom,isort,n,ind,j,i,k0,k,ig,nei,jatom,jsort,n1,ir,
     &           jnd
      real*8 :: phase,pi2
      complex*16 :: cc,c1
      complex*16, allocatable :: pb0(:,:)
      allocate(pb0(n_pbmt,n_pbmt))
      pi2=pi+pi
      do k0=1,k_star(iq)
        k=k_list(k0,iq)
        ig=k_group(k)
        pb0=pb
        call sym_w_mt(ig,pnt(1,k),pb0,n_pbmt)
        do isort=1,nsort
          iatom=iat_1(isort)
          n=n_pbmt0(isort)-1
          ind=iopb(iatom)
          do nei=1,nn_sites(iatom)
            jatom=atom_nn(nei,iatom)
            jsort=is(jatom)
            ir=lat_nn(nei,iatom)
            n1=n_pbmt0(jsort)-1
            jnd=iopb(jatom)
            phase=pi2*dot_product(pnt(:,k),rbs(:,ir))
            c1=dcmplx(cos(phase),sin(phase))/nqdiv
            cc=(0.d0,0.d0)
            do j=jnd,jnd+n1
              do i=ind,ind+n
                cc=cc+mtp(j)*pb(j,i)*mtp(i)
              enddo
            enddo
            chi_neig(nei,isort)=chi_neig(nei,isort)+cc*c1
          enddo
        enddo
      enddo
      deallocate(pb0)
      end
